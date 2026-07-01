#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <vector>

namespace cipc {

// Axis-aligned bounding box.
struct AABB {
    Eigen::Vector3d min;
    Eigen::Vector3d max;

    AABB()
        : min(1e32, 1e32, 1e32)
        , max(-1e32, -1e32, -1e32)
    {
    }

    AABB(const Eigen::Vector3d& mn, const Eigen::Vector3d& mx)
        : min(mn)
        , max(mx)
    {
    }

    bool isEmpty() const
    {
        return (min.array() > max.array()).any();
    }

    void setEmpty()
    {
        min.setConstant(1e32);
        max.setConstant(-1e32);
    }

    void extend(const Eigen::Vector3d& p)
    {
        min = min.cwiseMin(p);
        max = max.cwiseMax(p);
    }

    void extend(const AABB& other)
    {
        min = min.cwiseMin(other.min);
        max = max.cwiseMax(other.max);
    }

    void expand(double eps)
    {
        min.array() -= eps;
        max.array() += eps;
    }

    bool overlaps(const AABB& other) const
    {
        return (min.array() <= other.max.array()).all() && (max.array() >= other.min.array()).all();
    }

    bool contains(const Eigen::Vector3d& p) const
    {
        return (p.array() >= min.array()).all() && (p.array() <= max.array()).all();
    }

    double surfaceArea() const
    {
        const Eigen::Vector3d d = max - min;
        return 2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0]);
    }

    Eigen::Vector3d center() const { return 0.5 * (min + max); }
    Eigen::Vector3d extent() const { return max - min; }
};

// Lightweight bounding-volume hierarchy over indexed primitives.
// The caller provides per-primitive AABBs; the BVH stores primitive indices and
// a hierarchy of internal-node/leaf-node bounding boxes.
class BVH {
public:
    struct Node {
        AABB bbox;
        int left = -1; // left child index for internal nodes, -1 for leaf
        int right = -1; // right child index for internal nodes
        int begin = -1; // primitive range [begin, end) for leaf
        int end = -1;

        bool isLeaf() const { return left < 0; }
    };

    BVH() = default;

    // Build from per-primitive AABBs. If primIds is empty, uses 0, 1, ..., n-1.
    void build(const std::vector<AABB>& primBBoxes, const std::vector<int>& primIds = {})
    {
        nodes_.clear();
        primIds_.clear();

        const int n = static_cast<int>(primBBoxes.size());
        if (n == 0) {
            return;
        }

        std::vector<int> ids;
        if (primIds.empty()) {
            ids.resize(n);
            for (int i = 0; i < n; ++i) {
                ids[i] = i;
            }
        } else {
            ids = primIds;
        }

        std::vector<BuildItem> items;
        items.reserve(ids.size());
        for (int pid : ids) {
            items.push_back({ pid, primBBoxes[pid], primBBoxes[pid].center() });
        }

        buildNode(items, 0, static_cast<int>(items.size()));

        primIds_.resize(items.size());
        for (size_t i = 0; i < items.size(); ++i) {
            primIds_[i] = items[i].primId;
        }
    }

    // Refit leaf and internal node AABBs using new per-primitive AABBs.
    // The primitive identity and leaf membership must match the order used in build().
    void refit(const std::vector<AABB>& primBBoxes)
    {
        if (nodes_.empty()) {
            return;
        }
        refitNode(0, primBBoxes);
    }

    int rootIndex() const { return nodes_.empty() ? -1 : 0; }
    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<int>& primitiveIds() const { return primIds_; }

    // Query all primitives whose AABB overlaps queryBox. Callback receives primitive id.
    template <typename Callback>
    void queryAABB(const AABB& queryBox, Callback&& cb) const
    {
        if (nodes_.empty()) {
            return;
        }

        std::array<int, 256> stack;
        int stackPos = 0;
        stack[stackPos++] = 0;

        while (stackPos > 0) {
            const int nodeIdx = stack[--stackPos];
            const Node& node = nodes_[nodeIdx];
            if (!node.bbox.overlaps(queryBox)) {
                continue;
            }

            if (node.isLeaf()) {
                for (int i = node.begin; i < node.end; ++i) {
                    cb(primIds_[i]);
                }
            } else {
                if (node.left >= 0) {
                    stack[stackPos++] = node.left;
                }
                if (node.right >= 0) {
                    stack[stackPos++] = node.right;
                }
            }
        }
    }

    // Query all primitives whose AABB overlaps the swept AABB from aabb0 to aabb1.
    template <typename Callback>
    void querySweptAABB(const AABB& aabb0, const AABB& aabb1, Callback&& cb) const
    {
        AABB queryBox = aabb0;
        queryBox.extend(aabb1);
        queryAABB(queryBox, cb);
    }

private:
    struct BuildItem {
        int primId;
        AABB bbox;
        Eigen::Vector3d centroid;
    };

    static constexpr int kLeafCapacity = 4;
    static constexpr int kMaxDepth = 64;

    std::vector<Node> nodes_;
    std::vector<int> primIds_; // primitive ids reordered into leaf order

    int buildNode(std::vector<BuildItem>& items, int begin, int end, int depth = 0)
    {
        const int nodeIdx = static_cast<int>(nodes_.size());
        nodes_.emplace_back();

        AABB nodeBBox;
        AABB centroidBBox;
        for (int i = begin; i < end; ++i) {
            nodeBBox.extend(items[i].bbox);
            centroidBBox.extend(items[i].centroid);
        }
        nodes_[nodeIdx].bbox = nodeBBox;

        const int count = end - begin;
        if (count <= kLeafCapacity || depth >= kMaxDepth) {
            nodes_[nodeIdx].begin = begin;
            nodes_[nodeIdx].end = end;
            return nodeIdx;
        }

        // Choose split axis as the largest extent of the centroid bbox.
        const Eigen::Vector3d centroidExtent = centroidBBox.extent();
        int axis = 0;
        if (centroidExtent[1] > centroidExtent[axis]) {
            axis = 1;
        }
        if (centroidExtent[2] > centroidExtent[axis]) {
            axis = 2;
        }

        // Split at the median of the centroids along the chosen axis.
        const int mid = begin + count / 2;
        std::nth_element(items.begin() + begin, items.begin() + mid, items.begin() + end,
            [axis](const BuildItem& a, const BuildItem& b) {
                return a.centroid[axis] < b.centroid[axis];
            });

        const int left = buildNode(items, begin, mid, depth + 1);
        const int right = buildNode(items, mid, end, depth + 1);
        nodes_[nodeIdx].left = left;
        nodes_[nodeIdx].right = right;
        return nodeIdx;
    }

    void refitNode(int nodeId, const std::vector<AABB>& primBBoxes)
    {
        Node& node = nodes_[nodeId];
        if (node.isLeaf()) {
            node.bbox.setEmpty();
            for (int i = node.begin; i < node.end; ++i) {
                node.bbox.extend(primBBoxes[primIds_[i]]);
            }
            return;
        }

        if (node.left >= 0) {
            refitNode(node.left, primBBoxes);
        }
        if (node.right >= 0) {
            refitNode(node.right, primBBoxes);
        }

        node.bbox.setEmpty();
        if (node.left >= 0) {
            node.bbox.extend(nodes_[node.left].bbox);
        }
        if (node.right >= 0) {
            node.bbox.extend(nodes_[node.right].bbox);
        }
    }
};

} // namespace cipc
