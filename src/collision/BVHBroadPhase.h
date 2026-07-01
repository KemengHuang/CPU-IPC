#pragma once

#include "collision/BVH.h"
#include "mesh/Mesh.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>

// BVH-based replacement for SpatialHash. The public interface mirrors
// SpatialHash so that existing callers can switch via a compile-time alias.
class BVHBroadPhase {
public:
    BVHBroadPhase() = default;
    BVHBroadPhase(const mesh3D& mesh, double voxelSize, bool use_V_prev = false);
    BVHBroadPhase(const mesh3D& mesh, const Eigen::VectorXd& searchDir,
        double curMaxStepSize, double voxelSize, bool use_V_prev = false);

    void build(const mesh3D& mesh, double voxelSize);
    void build(const mesh3D& mesh, const std::vector<Eigen::Vector3d>& searchDir,
        double& curMaxStepSize, double voxelSize, bool use_V_prev = false);

    // Voxel-specific helpers are no-ops for BVH (kept for interface compatibility).
    void locateVoxelAxisIndex(const Eigen::Vector3d& pos, Eigen::Array<int, 1, 3>& voxelAxisIndex) const;
    int locateVoxelIndex(const Eigen::Vector3d& pos) const;
    int voxelAxisIndex2VoxelIndex(const int voxelAxisIndex[3]) const;
    int voxelAxisIndex2VoxelIndex(int ix, int iy, int iz) const;

    void queryPointForTriangles(const Eigen::Vector3d& pos, double radius, std::unordered_set<int>& triInds) const;
    void queryPointForTriangles(int svI, std::unordered_set<int>& sTriInds) const;
    void queryTriangleForPoints(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        double radius, std::unordered_set<int>& pointInds) const;
    void queryTriangleForEdges(const Eigen::Vector3d& v0,
        const Eigen::Vector3d& v1,
        const Eigen::Vector3d& v2,
        double radius, std::unordered_set<int>& edgeInds) const;
    void queryPointForPrimitives(const Eigen::Vector3d& pos,
        const Eigen::Vector3d& dir, std::unordered_set<int>& sVInds,
        std::unordered_set<int>& sEdgeInds, std::unordered_set<int>& sTriInds) const;
    void queryEdgeForPE(const Eigen::Vector3d& vBegin,
        const Eigen::Vector3d& vEnd,
        std::vector<int>& svInds, std::vector<int>& edgeInds) const;
    void queryEdgeForEdgesWithBBoxCheck(
        const mesh3D& mesh,
        const Eigen::Matrix<double, 1, 3>& vBegin,
        const Eigen::Matrix<double, 1, 3>& vEnd,
        double radius, std::vector<int>& edgeInds,
        int eIq = -1) const;
    void queryEdgeForEdgesWithBBoxCheck(const mesh3D& mesh,
        const std::vector<Eigen::Vector3d>& searchDir, double curMaxStepSize,
        int seI, std::unordered_set<int>& sEdgeInds) const;

    void calculateActivateSet(mesh3D& mesh);

private:
    cipc::BVH triBVH_; // over surface triangles
    cipc::BVH edgeBVH_; // over surface edges

    // Per-primitive AABBs stored in the *primitive index* order (not leaf order).
    std::vector<cipc::AABB> triAABBs_;
    std::vector<cipc::AABB> edgeAABBs_;
    std::vector<cipc::AABB> pointAABBs_;

    // For motion-aware queries we need to know which vertex positions were used.
    bool isMotionBuild_ = false;
    double motionStepSize_ = 0.0;
    std::vector<Eigen::Vector3d> motionSearchDir_;

    void collectPointTrianglePairs(mesh3D& mesh, std::vector<MMCVID>& ptPairs);
    void collectEdgeEdgePairs(mesh3D& mesh, std::vector<std::vector<MMCVID>>& eePairsPerEdge);
    void deduplicatePPPEPairs(mesh3D& mesh,
        const std::vector<MMCVID>& ptPairs,
        const std::vector<std::vector<MMCVID>>& eePairsPerEdge);

    static cipc::AABB triangleAABB(const mesh3D& mesh, int sfI);
    static cipc::AABB triangleSweptAABB(const mesh3D& mesh, int sfI,
        const std::vector<Eigen::Vector3d>& V,
        const std::vector<Eigen::Vector3d>& searchDir, double stepSize);
    static cipc::AABB edgeAABB(const mesh3D& mesh, int seI);
    static cipc::AABB edgeSweptAABB(const mesh3D& mesh, int seI,
        const std::vector<Eigen::Vector3d>& V,
        const std::vector<Eigen::Vector3d>& searchDir, double stepSize);
};
