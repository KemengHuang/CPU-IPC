#pragma once

#include "SimulationMesh.h"

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

struct CPUAABB {
    Eigen::Vector3d lower;
    Eigen::Vector3d upper;

    CPUAABB();
    void combine(const Eigen::Vector3d& point);
    void combine(const CPUAABB& box);
    Eigen::Vector3d center() const;
    bool overlaps(const CPUAABB& box, double gap = 0.0) const;
};

class CPULinearBVH {
public:
    void build(const std::vector<CPUAABB>& primitiveBoxes);
    void query(const CPUAABB& box, double gap, std::vector<int>& primitiveIds) const;

    bool empty() const { return nodes_.empty(); }
    std::size_t primitiveCount() const { return primitiveCount_; }

private:
    static constexpr std::uint32_t invalidIndex = 0xFFFFFFFFu;

    struct Node {
        std::uint32_t parent = invalidIndex;
        std::uint32_t left = invalidIndex;
        std::uint32_t right = invalidIndex;
        std::uint32_t primitive = invalidIndex;

        bool isLeaf() const { return primitive != invalidIndex; }
    };

    std::vector<CPUAABB> nodeBoxes_;
    std::vector<Node> nodes_;
    std::vector<std::pair<std::uint64_t, std::uint32_t>> mortonEntries_;
    std::vector<std::uint64_t> sortedKeys_;
    std::vector<std::pair<std::uint32_t, bool>> buildStack_;
    std::size_t primitiveCount_ = 0;
};

// CPU port of GPU_IPC/GPU_IPC/mlbvh.cu: separate face/edge LBVHs,
// centroid Morton keys, Karras topology, and stack-based PT/EE traversal.
class CPUBroadPhaseBVH {
public:
    void build(const mesh3D& mesh);
    void buildSwept(
        const mesh3D& mesh,
        const std::vector<Eigen::Vector3d>& searchDirection,
        double stepSize);

    void queryPointForTriangles(
        const Eigen::Vector3d& position,
        double radius,
        std::vector<int>& triangleIds) const;
    void querySweptPointForTriangles(
        int surfaceVertexIndex,
        double radius,
        std::vector<int>& triangleIds) const;
    void queryEdgeForEdges(
        const CPUAABB& edgeBox,
        double radius,
        std::vector<int>& edgeIds) const;
    void querySweptEdgeForEdges(
        int edgeIndex,
        double radius,
        std::vector<int>& edgeIds) const;

private:
    CPULinearBVH faceTree_;
    CPULinearBVH edgeTree_;
    std::vector<CPUAABB> sweptPointBoxes_;
    std::vector<CPUAABB> faceBoxes_;
    std::vector<CPUAABB> edgeBoxes_;
};

CPUAABB makePointAABB(const Eigen::Vector3d& point);
CPUAABB makeEdgeAABB(const Eigen::Vector3d& first, const Eigen::Vector3d& second);
CPUAABB makeTriangleAABB(
    const Eigen::Vector3d& first,
    const Eigen::Vector3d& second,
    const Eigen::Vector3d& third);
CPUAABB makeSweptPointAABB(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& searchDirection,
    double stepSize);
CPUAABB makeSweptEdgeAABB(
    const Eigen::Vector3d& first,
    const Eigen::Vector3d& firstDirection,
    const Eigen::Vector3d& second,
    const Eigen::Vector3d& secondDirection,
    double stepSize);
CPUAABB makeSweptTriangleAABB(
    const std::array<Eigen::Vector3d, 3>& positions,
    const std::array<Eigen::Vector3d, 3>& searchDirections,
    double stepSize);
