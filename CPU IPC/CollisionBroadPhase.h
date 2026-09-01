#pragma once
#ifndef CIPC_COLLISION_BROAD_PHASE_H
#define CIPC_COLLISION_BROAD_PHASE_H

#include "SimulationMesh.h"
#include "LBVH.h"
#include <unordered_map>
#include <unordered_set>
using namespace Eigen;

enum class BroadPhaseBackend {
    SpatialHash,
    LinearBVH
};

class SpatialHash {
public: // data
    Vector3d leftBottomCorner, rightTopCorner;
    double one_div_voxelSize;
    Array<int, 1, 3> voxelCount;
    int voxelCount0x1;

    int surfEdgeStartInd, surfTriStartInd;

    std::unordered_map<int, std::vector<int>> voxel;
    std::vector<std::vector<int>> pointAndEdgeOccupancy;
    std::vector<Eigen::Array<int, 1, 3>> surfaceVoxelAxisIndices;
    std::vector<Eigen::Array<int, 1, 3>> sweptMinVoxelAxisIndices;
    std::vector<Eigen::Array<int, 1, 3>> sweptMaxVoxelAxisIndices;
    std::vector<int> vertexToSurfaceIndex;
    std::vector<Vector3d> sweptSurfacePositions;
    std::vector<std::vector<int>> edgeVoxelLocations;
    std::vector<std::vector<int>> faceVoxelLocations;
    std::vector<int> activeVoxelKeys;
    BroadPhaseBackend backend = BroadPhaseBackend::SpatialHash;
    CPUBroadPhaseBVH bvh;
    double sweptQueryRadius = 0.0;

public: // constructor
    SpatialHash() = default;

    void setBackend(BroadPhaseBackend selectedBackend) { backend = selectedBackend; }
    BroadPhaseBackend getBackend() const { return backend; }

    void build(const mesh3D& mesh, double voxelSize);
    void build(const  mesh3D& mesh, const vector<Vector3d>& searchDir, double& curMaxStepSize, double voxelSize, bool use_V_prev = false);

    void locateVoxelAxisIndex(const Vector3d& pos, Eigen::Array<int, 1, 3>& voxelAxisIndex);
    int locateVoxelIndex(const Vector3d& pos);
    int voxelAxisIndex2VoxelIndex(const int voxelAxisIndex[3]);
    int voxelAxisIndex2VoxelIndex(int ix, int iy, int iz);
    void resetVoxelStorage();
    void insertVoxelEntry(int voxelIndex, int primitiveIndex);
    void queryPointForTriangles(const Vector3d& pos, double radius, std::vector<int>& triInds);
    void queryPointForTriangles(int svI, std::vector<int>& sTriInds);
    void queryTriangleForPoints(const Vector3d& v0,
        const Vector3d& v1,
        const Vector3d& v2,
        double radius, std::unordered_set<int>& pointInds);
    void queryTriangleForEdges(const Vector3d& v0,
        const Vector3d& v1,
        const Vector3d& v2,
        double radius, std::vector<int>& edgeInds);
    void queryPointForPrimitives(const Vector3d& pos,
        const Vector3d& dir, std::unordered_set<int>& sVInds,
        std::unordered_set<int>& sEdgeInds, std::unordered_set<int>& sTriInds);
    void queryEdgeForPE(const Vector3d& vBegin,
        const Vector3d& vEnd,
        std::vector<int>& svInds, std::vector<int>& edgeInds);
    void queryEdgeForEdgesWithBBoxCheck(
        const mesh3D& mesh,
        const Eigen::Matrix<double, 1, 3>& vBegin,
        const Eigen::Matrix<double, 1, 3>& vEnd,
        double radius, std::vector<int>& edgeInds,
        int eIq = -1);
    void queryEdgeForEdgesWithBBoxCheck(const mesh3D& mesh,
        const vector<Eigen::Vector3d>& searchDir, double curMaxStepSize,
        int seI, std::vector<int>& sEdgeInds);
    void calculateActivateSet(mesh3D& mesh);
};

class Ground {
public: // data
    Vector3d normal;
    double D;
public:
    Ground();
    void init(const Vector3d& m_normal, const double& d);
    double calculateGapFromObj(const mesh3D& mesh, const int& vId) const;
    void calculateActivateSet(mesh3D& mesh) const;
    //void calculateConstraintVal();
};




#endif // CIPC_COLLISION_BROAD_PHASE_H

