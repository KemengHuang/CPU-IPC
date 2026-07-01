#pragma once

#include "mesh/Mesh.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>

class SpatialHash {
public: // data
    Eigen::Vector3d leftBottomCorner, rightTopCorner;
    double one_div_voxelSize;
    Eigen::Array<int, 1, 3> voxelCount;
    int voxelCount0x1;

    int surfEdgeStartInd, surfTriStartInd;

    std::unordered_map<int, std::vector<int>> voxel;
    std::vector<std::vector<int>> pointAndEdgeOccupancy;

public: // constructor
    SpatialHash(void) {}
    SpatialHash(const mesh3D& mesh, double voxelSize, bool use_V_prev = false);
    SpatialHash(const mesh3D& mesh, const Eigen::VectorXd& searchDir, double curMaxStepSize, double voxelSize, bool use_V_prev = false);

    void build(const mesh3D& mesh, double voxelSize);
    void build(const mesh3D& mesh, const std::vector<Eigen::Vector3d>& searchDir, double& curMaxStepSize, double voxelSize, bool use_V_prev = false);

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
    void computeVoxelGrid(const mesh3D& mesh, double voxelSize);
    void collectPointTrianglePairs(mesh3D& mesh, std::vector<MMCVID>& ptPairs);
    void collectEdgeEdgePairs(mesh3D& mesh, std::vector<std::vector<MMCVID>>& eePairsPerEdge);
    void deduplicatePPPEPairs(mesh3D& mesh,
        const std::vector<MMCVID>& ptPairs,
        const std::vector<std::vector<MMCVID>>& eePairsPerEdge);
};

class Ground {
public: // data
    Eigen::Vector3d normal;
    double D;
public:
    Ground();
    void init(const Eigen::Vector3d& m_normal, const double& d);
    double calculateGapFromObj(const mesh3D& mesh, const int& vId) const;
    void calculateActivateSet(mesh3D& mesh) const;
};
