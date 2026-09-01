#include "LBVH.h"

// CPU/TBB adaptation of the Morton/Karras topology used by
// GPU_IPC/GPU_IPC/mlbvh.cu (Kemeng Huang, 2022-2024).

#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

std::uint32_t expandBits(std::uint32_t value)
{
    value = (value * 0x00010001u) & 0xFF0000FFu;
    value = (value * 0x00000101u) & 0x0F00F00Fu;
    value = (value * 0x00000011u) & 0xC30C30C3u;
    value = (value * 0x00000005u) & 0x49249249u;
    return value;
}

std::uint32_t mortonCode(double x, double y, double z)
{
    constexpr double resolution = 1024.0;
    auto quantize = [=](double value) {
        value = (std::max)(0.0, (std::min)(value * resolution, resolution - 1.0));
        return static_cast<std::uint32_t>(value);
    };
    const std::uint32_t xx = expandBits(quantize(x));
    const std::uint32_t yy = expandBits(quantize(y));
    const std::uint32_t zz = expandBits(quantize(z));
    return (xx << 2) | (yy << 1) | zz;
}

int countLeadingZeros(std::uint64_t value)
{
    if (value == 0) {
        return 64;
    }
    int result = 0;
    for (std::uint64_t mask = std::uint64_t(1) << 63; (value & mask) == 0; mask >>= 1) {
        ++result;
    }
    return result;
}

int commonPrefix(
    const std::vector<std::uint64_t>& keys,
    int first,
    int second)
{
    if (second < 0 || second >= static_cast<int>(keys.size())) {
        return -1;
    }
    return countLeadingZeros(keys[first] ^ keys[second]);
}

std::pair<std::uint32_t, std::uint32_t> determineRange(
    const std::vector<std::uint64_t>& keys,
    std::uint32_t index)
{
    const int leafCount = static_cast<int>(keys.size());
    if (index == 0) {
        return { 0u, static_cast<std::uint32_t>(leafCount - 1) };
    }

    const int current = static_cast<int>(index);
    const int leftPrefix = commonPrefix(keys, current, current - 1);
    const int rightPrefix = commonPrefix(keys, current, current + 1);
    const int direction = rightPrefix > leftPrefix ? 1 : -1;
    const int minimumPrefix = (std::min)(leftPrefix, rightPrefix);

    int maximumLength = 2;
    while (commonPrefix(keys, current, current + direction * maximumLength) > minimumPrefix) {
        maximumLength <<= 1;
    }

    int length = 0;
    for (int step = maximumLength >> 1; step > 0; step >>= 1) {
        if (commonPrefix(keys, current, current + direction * (length + step)) > minimumPrefix) {
            length += step;
        }
    }

    const int other = current + direction * length;
    return {
        static_cast<std::uint32_t>((std::min)(current, other)),
        static_cast<std::uint32_t>((std::max)(current, other))
    };
}

std::uint32_t findSplit(
    const std::vector<std::uint64_t>& keys,
    std::uint32_t first,
    std::uint32_t last)
{
    const std::uint64_t firstCode = keys[first];
    const std::uint64_t lastCode = keys[last];
    if (firstCode == lastCode) {
        return (first + last) >> 1;
    }

    const int common = countLeadingZeros(firstCode ^ lastCode);
    std::uint32_t split = first;
    std::uint32_t step = last - first;
    do {
        step = (step + 1) >> 1;
        const std::uint32_t candidate = split + step;
        if (candidate < last
            && countLeadingZeros(firstCode ^ keys[candidate]) > common) {
            split = candidate;
        }
    } while (step > 1);
    return split;
}

CPUAABB sweptPointBox(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& searchDirection,
    double stepSize)
{
    CPUAABB box;
    box.combine(position);
    box.combine(position - stepSize * searchDirection);
    return box;
}

} // namespace

CPUAABB::CPUAABB()
    : lower(Eigen::Vector3d::Constant(std::numeric_limits<double>::max())),
      upper(Eigen::Vector3d::Constant(std::numeric_limits<double>::lowest()))
{
}

void CPUAABB::combine(const Eigen::Vector3d& point)
{
    lower = lower.cwiseMin(point);
    upper = upper.cwiseMax(point);
}

void CPUAABB::combine(const CPUAABB& box)
{
    lower = lower.cwiseMin(box.lower);
    upper = upper.cwiseMax(box.upper);
}

Eigen::Vector3d CPUAABB::center() const
{
    return 0.5 * (lower + upper);
}

bool CPUAABB::overlaps(const CPUAABB& box, double gap) const
{
    return !((box.lower.array() - upper.array() > gap).any()
        || (lower.array() - box.upper.array() > gap).any());
}

void CPULinearBVH::build(const std::vector<CPUAABB>& primitiveBoxes)
{
    primitiveCount_ = primitiveBoxes.size();
    const std::size_t leafCount = primitiveBoxes.size();
    if (leafCount > std::numeric_limits<std::uint32_t>::max()) {
        throw std::overflow_error("CPU LBVH primitive count exceeds uint32 range");
    }
    if (leafCount == 0) {
        nodes_.clear();
        nodeBoxes_.clear();
        return;
    }

    CPUAABB scene;
    for (const CPUAABB& box : primitiveBoxes) {
        scene.combine(box);
    }
    const Eigen::Vector3d extent = scene.upper - scene.lower;

    mortonEntries_.resize(leafCount);
    tbb::parallel_for(std::size_t(0), leafCount, [&](std::size_t primitive) {
        const Eigen::Vector3d offset = primitiveBoxes[primitive].center() - scene.lower;
        const double x = extent.x() > 0.0 ? offset.x() / extent.x() : 0.0;
        const double y = extent.y() > 0.0 ? offset.y() / extent.y() : 0.0;
        const double z = extent.z() > 0.0 ? offset.z() / extent.z() : 0.0;
        const std::uint32_t code = mortonCode(x, y, z);
        mortonEntries_[primitive] = {
            (static_cast<std::uint64_t>(code) << 32) | primitive,
            static_cast<std::uint32_t>(primitive)
        };
    });
    tbb::parallel_sort(mortonEntries_.begin(), mortonEntries_.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });

    const std::size_t nodeCount = 2 * leafCount - 1;
    nodes_.assign(nodeCount, Node{});
    nodeBoxes_.assign(nodeCount, CPUAABB{});
    sortedKeys_.resize(leafCount);
    for (std::size_t sorted = 0; sorted < leafCount; ++sorted) {
        sortedKeys_[sorted] = mortonEntries_[sorted].first;
        const std::size_t leaf = sorted + leafCount - 1;
        nodes_[leaf].primitive = mortonEntries_[sorted].second;
        nodeBoxes_[leaf] = primitiveBoxes[mortonEntries_[sorted].second];
    }

    if (leafCount == 1) {
        return;
    }

    tbb::parallel_for(std::size_t(0), leafCount - 1, [&](std::size_t internal) {
        const auto range = determineRange(sortedKeys_, static_cast<std::uint32_t>(internal));
        const std::uint32_t split = findSplit(sortedKeys_, range.first, range.second);
        std::uint32_t left = split;
        std::uint32_t right = split + 1;
        if (split == range.first) {
            left += static_cast<std::uint32_t>(leafCount - 1);
        }
        if (split + 1 == range.second) {
            right += static_cast<std::uint32_t>(leafCount - 1);
        }
        nodes_[internal].left = left;
        nodes_[internal].right = right;
        nodes_[left].parent = static_cast<std::uint32_t>(internal);
        nodes_[right].parent = static_cast<std::uint32_t>(internal);
    });

    buildStack_.clear();
    buildStack_.reserve(nodeCount);
    buildStack_.emplace_back(0u, false);
    while (!buildStack_.empty()) {
        const auto current = buildStack_.back();
        buildStack_.pop_back();
        const Node& node = nodes_[current.first];
        if (node.isLeaf()) {
            continue;
        }
        if (current.second) {
            nodeBoxes_[current.first] = nodeBoxes_[node.left];
            nodeBoxes_[current.first].combine(nodeBoxes_[node.right]);
        }
        else {
            buildStack_.emplace_back(current.first, true);
            buildStack_.emplace_back(node.right, false);
            buildStack_.emplace_back(node.left, false);
        }
    }
}

void CPULinearBVH::query(
    const CPUAABB& box,
    double gap,
    std::vector<int>& primitiveIds) const
{
    primitiveIds.clear();
    if (nodes_.empty()) {
        return;
    }

    std::array<std::uint32_t, 128> stack;
    std::size_t size = 0;
    stack[size++] = 0;
    while (size > 0) {
        const std::uint32_t nodeIndex = stack[--size];
        if (!box.overlaps(nodeBoxes_[nodeIndex], gap)) {
            continue;
        }
        const Node& node = nodes_[nodeIndex];
        if (node.isLeaf()) {
            primitiveIds.emplace_back(static_cast<int>(node.primitive));
            continue;
        }
        if (size + 2 > stack.size()) {
            throw std::runtime_error("CPU LBVH traversal stack overflow");
        }
        stack[size++] = node.right;
        stack[size++] = node.left;
    }
    std::sort(primitiveIds.begin(), primitiveIds.end());
}

void CPUBroadPhaseBVH::build(const mesh3D& mesh)
{
    faceBoxes_.resize(mesh.surface.size());
    tbb::parallel_for(std::size_t(0), mesh.surface.size(), [&](std::size_t face) {
        const Vector4i& vertices = mesh.surface[face];
        faceBoxes_[face] = makeTriangleAABB(
            mesh.vertexes[vertices[0]],
            mesh.vertexes[vertices[1]],
            mesh.vertexes[vertices[2]]);
    });

    edgeBoxes_.resize(mesh.surfEdges.size());
    tbb::parallel_for(std::size_t(0), mesh.surfEdges.size(), [&](std::size_t edge) {
        const auto& vertices = mesh.surfEdges[edge];
        edgeBoxes_[edge] = makeEdgeAABB(
            mesh.vertexes[vertices.first],
            mesh.vertexes[vertices.second]);
    });
    sweptPointBoxes_.clear();
    faceTree_.build(faceBoxes_);
    edgeTree_.build(edgeBoxes_);
}

void CPUBroadPhaseBVH::buildSwept(
    const mesh3D& mesh,
    const std::vector<Eigen::Vector3d>& searchDirection,
    double stepSize)
{
    sweptPointBoxes_.resize(mesh.surfVerts.size());
    tbb::parallel_for(std::size_t(0), mesh.surfVerts.size(), [&](std::size_t surfaceVertex) {
        const std::size_t vertex = static_cast<std::size_t>(mesh.surfVerts[surfaceVertex]);
        sweptPointBoxes_[surfaceVertex] = makeSweptPointAABB(
            mesh.vertexes[vertex], searchDirection[vertex], stepSize);
    });

    faceBoxes_.resize(mesh.surface.size());
    tbb::parallel_for(std::size_t(0), mesh.surface.size(), [&](std::size_t face) {
        const Vector4i& vertices = mesh.surface[face];
        faceBoxes_[face] = makeSweptTriangleAABB(
            { mesh.vertexes[vertices[0]], mesh.vertexes[vertices[1]], mesh.vertexes[vertices[2]] },
            { searchDirection[vertices[0]], searchDirection[vertices[1]], searchDirection[vertices[2]] },
            stepSize);
    });

    edgeBoxes_.resize(mesh.surfEdges.size());
    tbb::parallel_for(std::size_t(0), mesh.surfEdges.size(), [&](std::size_t edge) {
        const auto& vertices = mesh.surfEdges[edge];
        edgeBoxes_[edge] = makeSweptEdgeAABB(
            mesh.vertexes[vertices.first], searchDirection[vertices.first],
            mesh.vertexes[vertices.second], searchDirection[vertices.second],
            stepSize);
    });

    faceTree_.build(faceBoxes_);
    edgeTree_.build(edgeBoxes_);
}

void CPUBroadPhaseBVH::queryPointForTriangles(
    const Eigen::Vector3d& position,
    double radius,
    std::vector<int>& triangleIds) const
{
    faceTree_.query(makePointAABB(position), radius, triangleIds);
}

void CPUBroadPhaseBVH::querySweptPointForTriangles(
    int surfaceVertexIndex,
    double radius,
    std::vector<int>& triangleIds) const
{
    faceTree_.query(sweptPointBoxes_.at(static_cast<std::size_t>(surfaceVertexIndex)), radius, triangleIds);
}

void CPUBroadPhaseBVH::queryEdgeForEdges(
    const CPUAABB& edgeBox,
    double radius,
    std::vector<int>& edgeIds) const
{
    edgeTree_.query(edgeBox, radius, edgeIds);
}

void CPUBroadPhaseBVH::querySweptEdgeForEdges(
    int edgeIndex,
    double radius,
    std::vector<int>& edgeIds) const
{
    edgeTree_.query(edgeBoxes_.at(static_cast<std::size_t>(edgeIndex)), radius, edgeIds);
}

CPUAABB makePointAABB(const Eigen::Vector3d& point)
{
    CPUAABB box;
    box.combine(point);
    return box;
}

CPUAABB makeEdgeAABB(const Eigen::Vector3d& first, const Eigen::Vector3d& second)
{
    CPUAABB box;
    box.combine(first);
    box.combine(second);
    return box;
}

CPUAABB makeTriangleAABB(
    const Eigen::Vector3d& first,
    const Eigen::Vector3d& second,
    const Eigen::Vector3d& third)
{
    CPUAABB box;
    box.combine(first);
    box.combine(second);
    box.combine(third);
    return box;
}

CPUAABB makeSweptPointAABB(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& searchDirection,
    double stepSize)
{
    return sweptPointBox(position, searchDirection, stepSize);
}

CPUAABB makeSweptEdgeAABB(
    const Eigen::Vector3d& first,
    const Eigen::Vector3d& firstDirection,
    const Eigen::Vector3d& second,
    const Eigen::Vector3d& secondDirection,
    double stepSize)
{
    CPUAABB box = makeSweptPointAABB(first, firstDirection, stepSize);
    box.combine(makeSweptPointAABB(second, secondDirection, stepSize));
    return box;
}

CPUAABB makeSweptTriangleAABB(
    const std::array<Eigen::Vector3d, 3>& positions,
    const std::array<Eigen::Vector3d, 3>& searchDirections,
    double stepSize)
{
    CPUAABB box;
    for (int local = 0; local < 3; ++local) {
        box.combine(makeSweptPointAABB(
            positions[local], searchDirections[local], stepSize));
    }
    return box;
}
