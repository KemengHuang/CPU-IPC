#include "ContactMechanics.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>
#include <tbb/enumerable_thread_specific.h>
#include "CollisionBroadPhase.h"
#include <algorithm>
#include <assert.h>
#include <stdexcept>
void SpatialHash::locateVoxelAxisIndex(const Vector3d& pos, Eigen::Array<int, 1, 3>& voxelAxisIndex)
{
    voxelAxisIndex = ((pos - leftBottomCorner) * one_div_voxelSize).array().floor().template cast<int>();
}

int SpatialHash::locateVoxelIndex(const Vector3d& pos)
{
    Eigen::Array<int, 1, 3> voxelAxisIndex;
    locateVoxelAxisIndex(pos, voxelAxisIndex);
    return voxelAxisIndex2VoxelIndex(voxelAxisIndex.data());
}

int SpatialHash::voxelAxisIndex2VoxelIndex(const int voxelAxisIndex[3])
{
    return voxelAxisIndex2VoxelIndex(voxelAxisIndex[0], voxelAxisIndex[1], voxelAxisIndex[2]);
}

int SpatialHash::voxelAxisIndex2VoxelIndex(int ix, int iy, int iz)
{
    return ix + iy * voxelCount[0] + iz * voxelCount0x1;
}

void SpatialHash::resetVoxelStorage()
{
    if (voxel.size() > activeVoxelKeys.size() * 4 + 1024) {
        voxel.clear();
    }
    else {
        for (const int key : activeVoxelKeys) {
            const auto entry = voxel.find(key);
            if (entry != voxel.end()) {
                entry->second.clear();
            }
        }
    }
    activeVoxelKeys.clear();
}

void SpatialHash::insertVoxelEntry(int voxelIndex, int primitiveIndex)
{
    std::vector<int>& entries = voxel[voxelIndex];
    if (entries.empty()) {
        activeVoxelKeys.emplace_back(voxelIndex);
    }
    entries.emplace_back(primitiveIndex);
}

void SpatialHash::build(const  mesh3D& mesh, double voxelSize)
{
    surfEdgeStartInd = static_cast<int>(mesh.surfVerts.size());
    surfTriStartInd = surfEdgeStartInd + static_cast<int>(mesh.surfEdges.size());
    sweptQueryRadius = 0.0;
    if (backend == BroadPhaseBackend::LinearBVH) {
        bvh.build(mesh);
        return;
    }

    //double xmin = DBL_MAX, ymin = DBL_MAX, zmin = DBL_MAX;
    //double xmax = DBL_MIN, ymax = DBL_MIN, zmax = DBL_MIN;

    leftBottomCorner = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, (int)mesh.surfVerts.size()), Vector3d(1e32, 1e32, 1e32), [&](const tbb::blocked_range<int>& rg, Vector3d sceneSize) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                Vector3d pos = mesh.vertexes[mesh.surfVerts[i]];
                sceneSize[0] = min(pos[0], sceneSize[0]);
                sceneSize[1] = min(pos[1], sceneSize[1]);
                sceneSize[2] = min(pos[2], sceneSize[2]);
            }
            return sceneSize;
        },
        [&](Vector3d left, Vector3d right) {
            Vector3d output;
            output[0] = min(left[0], right[0]);
            output[1] = min(left[1], right[1]);
            output[2] = min(left[2], right[2]);
            return output;
        }
        );

    rightTopCorner = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, (int)mesh.surfVerts.size()), Vector3d(-1e32, -1e32, -1e32), [&](const tbb::blocked_range<int>& rg, Vector3d sceneSize) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                Vector3d pos = mesh.vertexes[mesh.surfVerts[i]];
                sceneSize[0] = max(pos[0], sceneSize[0]);
                sceneSize[1] = max(pos[1], sceneSize[1]);
                sceneSize[2] = max(pos[2], sceneSize[2]);
            }
            return sceneSize;
        },
        [&](Vector3d left, Vector3d right) {
            Vector3d output;
            output[0] = max(left[0], right[0]);
            output[1] = max(left[1], right[1]);
            output[2] = max(left[2], right[2]);
            return output;
        }
        );
    //for (int i = 0;i < mesh.surface.size();i++) {
    //    for (int j = 0;j < 3;j++) {
    //        Vector3d pos = mesh.vertexes[mesh.surface[i][j]];
    //        if (xmin > pos[0]) xmin = pos[0];
    //        if (ymin > pos[1]) ymin = pos[1];
    //        if (zmin > pos[2]) zmin = pos[2];
    //        if (xmax < pos[0]) xmax = pos[0];
    //        if (ymax < pos[1]) ymax = pos[1];
    //        if (zmax < pos[2]) zmax = pos[2];
    //    }
    //}
    //leftBottomCorner = Vector3d(xmin, ymin, zmin);
    //rightTopCorner = Vector3d(xmax, ymax, zmax);

    one_div_voxelSize = 1.0 / voxelSize;
    Eigen::Array<double, 1, 3> range = rightTopCorner - leftBottomCorner;
    voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
    if (voxelCount.minCoeff() <= 0) {
        // cast overflow due to huge search direction
        one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
        voxelCount.setOnes();
    }
    voxelCount0x1 = voxelCount[0] * voxelCount[1];

    surfEdgeStartInd = mesh.surfVerts.size();
    surfTriStartInd = surfEdgeStartInd + mesh.surfEdges.size();

    surfaceVoxelAxisIndices.resize(mesh.surfVerts.size());
    vertexToSurfaceIndex.assign(mesh.vertexNum, -1);

    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI)
        //for (int svI = 0; svI < mesh.surfVerts.size(); ++svI)
        {
            int vI = mesh.surfVerts[svI];
            locateVoxelAxisIndex(mesh.vertexes[vI], surfaceVoxelAxisIndices[svI]);
            vertexToSurfaceIndex[vI] = svI;
        }
    );


    resetVoxelStorage();


    for (int svI = 0; svI < mesh.surfVerts.size(); ++svI) {
        insertVoxelEntry(locateVoxelIndex(mesh.vertexes[mesh.surfVerts[svI]]), svI);
    }

    edgeVoxelLocations.resize(mesh.surfEdges.size());
    for (auto& locations : edgeVoxelLocations) {
        locations.clear();
    }

    tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int seCount)
        //for (int seCount = 0; seCount < mesh.surfEdges.size(); ++seCount)
        {
            const auto& seI = mesh.surfEdges[seCount];

            const Eigen::Array<int, 1, 3>& voxelAxisIndex_first = surfaceVoxelAxisIndices[vertexToSurfaceIndex[seI.first]];
            const Eigen::Array<int, 1, 3>& voxelAxisIndex_second = surfaceVoxelAxisIndices[vertexToSurfaceIndex[seI.second]];
            Eigen::Array<int, 1, 3> mins = voxelAxisIndex_first.min(voxelAxisIndex_second);
            Eigen::Array<int, 1, 3> maxs = voxelAxisIndex_first.max(voxelAxisIndex_second);
            for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                int zOffset = iz * voxelCount0x1;
                for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                    int yzOffset = iy * voxelCount[0] + zOffset;
                    for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                        edgeVoxelLocations[seCount].emplace_back(ix + yzOffset);
                    }
                }
            }
        }
    );

    faceVoxelLocations.resize(mesh.surface.size());
    for (auto& locations : faceVoxelLocations) {
        locations.clear();
    }

    tbb::parallel_for(0, (int)mesh.surface.size(), 1, [&](int sfI)
        //for (int sfI = 0; sfI < mesh.surface.size(); ++sfI)
        {
            const Eigen::Array<int, 1, 3>& voxelAxisIndex0 = surfaceVoxelAxisIndices[vertexToSurfaceIndex[mesh.surface[sfI][0]]];
            const Eigen::Array<int, 1, 3>& voxelAxisIndex1 = surfaceVoxelAxisIndices[vertexToSurfaceIndex[mesh.surface[sfI][1]]];
            const Eigen::Array<int, 1, 3>& voxelAxisIndex2 = surfaceVoxelAxisIndices[vertexToSurfaceIndex[mesh.surface[sfI][2]]];
            Eigen::Array<int, 1, 3> mins = voxelAxisIndex0.min(voxelAxisIndex1).min(voxelAxisIndex2);
            Eigen::Array<int, 1, 3> maxs = voxelAxisIndex0.max(voxelAxisIndex1).max(voxelAxisIndex2);
            for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                int zOffset = iz * voxelCount0x1;
                for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                    int yzOffset = iy * voxelCount[0] + zOffset;
                    for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                        faceVoxelLocations[sfI].emplace_back(ix + yzOffset);
                    }
                }
            }
        }
    );

    for (int seCount = 0; seCount < edgeVoxelLocations.size(); ++seCount) {
        for (const auto& voxelI : edgeVoxelLocations[seCount]) {
            insertVoxelEntry(voxelI, seCount + surfEdgeStartInd);
        }
    }
    for (int sfI = 0; sfI < faceVoxelLocations.size(); ++sfI) {
        for (const auto& voxelI : faceVoxelLocations[sfI]) {
            insertVoxelEntry(voxelI, sfI + surfTriStartInd);
        }
    }
}






void SpatialHash::build(const  mesh3D& mesh, const vector<Vector3d>& searchDir, double& curMaxStepSize, double voxelSize, bool use_V_prev)
{
    surfEdgeStartInd = static_cast<int>(mesh.surfVerts.size());
    surfTriStartInd = surfEdgeStartInd + static_cast<int>(mesh.surfEdges.size());
    // ACCD uses zero thickness; swept primitive AABBs are already conservative.
    // Keep the CPU solver's existing candidate semantics (the GPU code uses
    // sqrt(dHat) here, which is a safe but noticeably larger superset).
    sweptQueryRadius = 0.0;
    if (backend == BroadPhaseBackend::LinearBVH) {
        if (use_V_prev) {
            throw std::logic_error("CPU LBVH does not support use_V_prev swept builds");
        }
        bvh.buildSwept(mesh, searchDir, curMaxStepSize);
        return;
    }

    double pSize = 0;
    for (int svI = 0; svI < mesh.surfVerts.size(); ++svI) {
        int vI = mesh.surfVerts[svI];
        pSize += std::abs(searchDir[vI][0]);
        pSize += std::abs(searchDir[vI][1]);
        pSize += std::abs(searchDir[vI][2]);
    }
    pSize /= mesh.surfVerts.size() * 3;

    const double spanSize = curMaxStepSize * pSize / voxelSize;
    //if (spanSize > 1) {
    //    curMaxStepSize /= spanSize;
    //    // curMaxStepSize reduced for CCD spatial hash efficiency
    //}

    const vector<Vector3d>& V = use_V_prev ? mesh.V_prev : mesh.vertexes;
    sweptSurfacePositions.resize(mesh.surfVerts.size());
    vertexToSurfaceIndex.assign(mesh.vertexNum, -1);

    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI)
        //for (int svI = 0; svI < mesh.surfVerts.size(); ++svI) 
        {
            int vI = mesh.surfVerts[svI];
            sweptSurfacePositions[svI] = V[vI] - curMaxStepSize * searchDir[vI];
            vertexToSurfaceIndex[vI] = svI;
        }
    );


    leftBottomCorner = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, (int)mesh.surfVerts.size()), Vector3d(1e32, 1e32, 1e32), [&](const tbb::blocked_range<int>& rg, Vector3d sceneSize) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                Vector3d pos = V[mesh.surfVerts[i]];
                sceneSize[0] = min(pos[0], sceneSize[0]);
                sceneSize[1] = min(pos[1], sceneSize[1]);
                sceneSize[2] = min(pos[2], sceneSize[2]);

                pos = sweptSurfacePositions[i];
                sceneSize[0] = min(pos[0], sceneSize[0]);
                sceneSize[1] = min(pos[1], sceneSize[1]);
                sceneSize[2] = min(pos[2], sceneSize[2]);

            }
            return sceneSize;
        },
        [&](Vector3d left, Vector3d right) {
            Vector3d output;
            output[0] = min(left[0], right[0]);
            output[1] = min(left[1], right[1]);
            output[2] = min(left[2], right[2]);
            return output;
        }
        );

    rightTopCorner = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, (int)mesh.surfVerts.size()), Vector3d(-1e32, -1e32, -1e32), [&](const tbb::blocked_range<int>& rg, Vector3d sceneSize) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                Vector3d pos = V[mesh.surfVerts[i]];
                sceneSize[0] = max(pos[0], sceneSize[0]);
                sceneSize[1] = max(pos[1], sceneSize[1]);
                sceneSize[2] = max(pos[2], sceneSize[2]);

                pos = sweptSurfacePositions[i];
                sceneSize[0] = max(pos[0], sceneSize[0]);
                sceneSize[1] = max(pos[1], sceneSize[1]);
                sceneSize[2] = max(pos[2], sceneSize[2]);
            }
            return sceneSize;
        },
        [&](Vector3d left, Vector3d right) {
            Vector3d output;
            output[0] = max(left[0], right[0]);
            output[1] = max(left[1], right[1]);
            output[2] = max(left[2], right[2]);
            return output;
        }
        );

    //double xmin = std::numeric_limits<double>::max(), ymin = std::numeric_limits<double>::max(), zmin = std::numeric_limits<double>::max();
    //double xmax = DBL_MIN, ymax = DBL_MIN, zmax = DBL_MIN;
    //for (int i = 0;i < mesh.surface.size();i++) {
    //    for (int j = 0;j < 3;j++) {
    //        Vector3d pos = V[mesh.surface[i][j]];
    //        if (xmin > pos[0]) xmin = pos[0];
    //        if (ymin > pos[1]) ymin = pos[1];
    //        if (zmin > pos[2]) zmin = pos[2];
    //        if (xmax < pos[0]) xmax = pos[0];
    //        if (ymax < pos[1]) ymax = pos[1];
    //        if (zmax < pos[2]) zmax = pos[2];

    //        pos = SVt[vI2SVI[mesh.surface[i][j]]];
    //        if (xmin > pos[0]) xmin = pos[0];
    //        if (ymin > pos[1]) ymin = pos[1];
    //        if (zmin > pos[2]) zmin = pos[2];
    //        if (xmax < pos[0]) xmax = pos[0];
    //        if (ymax < pos[1]) ymax = pos[1];
    //        if (zmax < pos[2]) zmax = pos[2];
    //    }
    //}

    //leftBottomCorner = Vector3d(xmin, ymin, zmin);
    //rightTopCorner = Vector3d(xmax, ymax, zmax);

    one_div_voxelSize = 1.0 / voxelSize;
    Eigen::Array<double, 1, 3> range = rightTopCorner - leftBottomCorner;
    voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
    if (voxelCount.minCoeff() <= 0) {
        // cast overflow due to huge search direction
        one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
        voxelCount.setOnes();
    }
    voxelCount0x1 = voxelCount[0] * voxelCount[1];

    surfEdgeStartInd = mesh.surfVerts.size();
    surfTriStartInd = surfEdgeStartInd + mesh.surfEdges.size();

    // precompute svVAI
    sweptMinVoxelAxisIndices.resize(mesh.surfVerts.size());
    sweptMaxVoxelAxisIndices.resize(mesh.surfVerts.size());

    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI)
        {
            int vI = mesh.surfVerts[svI];
            Eigen::Array<int, 1, 3> v0VAI, vtVAI;
            locateVoxelAxisIndex(V[vI], v0VAI);
            locateVoxelAxisIndex(sweptSurfacePositions[svI], vtVAI);
            sweptMinVoxelAxisIndices[svI] = v0VAI.min(vtVAI);
            sweptMaxVoxelAxisIndices[svI] = v0VAI.max(vtVAI);
        }

    );

    resetVoxelStorage();
    pointAndEdgeOccupancy.resize(surfTriStartInd);
    for (auto& occupancy : pointAndEdgeOccupancy) {
        occupancy.clear();
    }


    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI)

        {
            const Eigen::Array<int, 1, 3>& mins = sweptMinVoxelAxisIndices[svI];
            const Eigen::Array<int, 1, 3>& maxs = sweptMaxVoxelAxisIndices[svI];
            pointAndEdgeOccupancy[svI].reserve((maxs - mins + 1).prod());
            for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                int zOffset = iz * voxelCount0x1;
                for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                    int yzOffset = iy * voxelCount[0] + zOffset;
                    for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                        pointAndEdgeOccupancy[svI].emplace_back(ix + yzOffset);
                    }
                }
            }
        }

    );



    tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int seCount)

        {
            int seIInd = seCount + surfEdgeStartInd;
            const auto& seI = mesh.surfEdges[seCount];

            Eigen::Array<int, 1, 3> mins = sweptMinVoxelAxisIndices[vertexToSurfaceIndex[seI.first]].min(sweptMinVoxelAxisIndices[vertexToSurfaceIndex[seI.second]]);
            Eigen::Array<int, 1, 3> maxs = sweptMaxVoxelAxisIndices[vertexToSurfaceIndex[seI.first]].max(sweptMaxVoxelAxisIndices[vertexToSurfaceIndex[seI.second]]);
            pointAndEdgeOccupancy[seIInd].reserve((maxs - mins + 1).prod());
            for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                int zOffset = iz * voxelCount0x1;
                for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                    int yzOffset = iy * voxelCount[0] + zOffset;
                    for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                        pointAndEdgeOccupancy[seIInd].emplace_back(ix + yzOffset);
                    }
                }
            }
        }

    );

    faceVoxelLocations.resize(mesh.surface.size());
    for (auto& locations : faceVoxelLocations) {
        locations.clear();
    }

    tbb::parallel_for(0, (int)mesh.surface.size(), 1, [&](int sfI)

        {
            Eigen::Array<int, 1, 3> mins = sweptMinVoxelAxisIndices[vertexToSurfaceIndex[mesh.surface[sfI][0]]]
                .min(sweptMinVoxelAxisIndices[vertexToSurfaceIndex[mesh.surface[sfI][1]]])
                .min(sweptMinVoxelAxisIndices[vertexToSurfaceIndex[mesh.surface[sfI][2]]]);
            Eigen::Array<int, 1, 3> maxs = sweptMaxVoxelAxisIndices[vertexToSurfaceIndex[mesh.surface[sfI][0]]]
                .max(sweptMaxVoxelAxisIndices[vertexToSurfaceIndex[mesh.surface[sfI][1]]])
                .max(sweptMaxVoxelAxisIndices[vertexToSurfaceIndex[mesh.surface[sfI][2]]]);
            for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                int zOffset = iz * voxelCount0x1;
                for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                    int yzOffset = iy * voxelCount[0] + zOffset;
                    for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                        faceVoxelLocations[sfI].emplace_back(ix + yzOffset);
                    }
                }
            }
        }

    );

    for (int i = 0; i < pointAndEdgeOccupancy.size(); ++i) {
        for (const auto& voxelI : pointAndEdgeOccupancy[i]) {
            insertVoxelEntry(voxelI, i);
        }
    }
    for (int sfI = 0; sfI < faceVoxelLocations.size(); ++sfI) {
        for (const auto& voxelI : faceVoxelLocations[sfI]) {
            insertVoxelEntry(voxelI, sfI + surfTriStartInd);
        }
    }
}




void SpatialHash::queryPointForTriangles(const Vector3d& pos, double radius, std::vector<int>& triInds)
{
    if (backend == BroadPhaseBackend::LinearBVH) {
        bvh.queryPointForTriangles(pos, radius, triInds);
        return;
    }
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(pos.array() - radius, mins);
    locateVoxelAxisIndex(pos.array() + radius, maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount.array() - 1);

    triInds.clear();
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= surfTriStartInd) {
                            triInds.emplace_back(indI - surfTriStartInd);
                        }
                    }
                }
            }
        }
    }
    std::sort(triInds.begin(), triInds.end());
    triInds.erase(std::unique(triInds.begin(), triInds.end()), triInds.end());
}

void SpatialHash::queryPointForTriangles(int svI, std::vector<int>& sTriInds)
{
    if (backend == BroadPhaseBackend::LinearBVH) {
        bvh.querySweptPointForTriangles(svI, sweptQueryRadius, sTriInds);
        return;
    }
    sTriInds.clear();
    for (const auto& voxelInd : pointAndEdgeOccupancy[svI]) {
        const auto& voxelI = voxel.find(voxelInd);
        assert(voxelI != voxel.end());
        for (const auto& indI : voxelI->second) {
            if (indI >= surfTriStartInd) {
                sTriInds.emplace_back(indI - surfTriStartInd);
            }
        }
    }
    std::sort(sTriInds.begin(), sTriInds.end());
    sTriInds.erase(std::unique(sTriInds.begin(), sTriInds.end()), sTriInds.end());
}

void SpatialHash::queryTriangleForPoints(const Vector3d& v0,
    const Vector3d& v1,
    const Vector3d& v2,
    double radius, std::unordered_set<int>& pointInds)
{
    Vector3d leftBottom = v0.array().min(v1.array()).min(v2.array());
    Vector3d rightTop = v0.array().max(v1.array()).max(v2.array());
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(leftBottom.array() - radius, mins);
    locateVoxelAxisIndex(rightTop.array() + radius, maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount - 1);

    pointInds.clear();
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI < surfEdgeStartInd) {
                            pointInds.insert(indI);
                        }
                    }
                }
            }
        }
    }
}



void SpatialHash::queryTriangleForEdges(const Vector3d& v0,
    const Vector3d& v1,
    const Vector3d& v2,
    double radius, std::vector<int>& edgeInds)
{
    if (backend == BroadPhaseBackend::LinearBVH) {
        bvh.queryEdgeForEdges(makeTriangleAABB(v0, v1, v2), radius, edgeInds);
        return;
    }
    Vector3d leftBottom = v0.array().min(v1.array()).min(v2.array());
    Vector3d rightTop = v0.array().max(v1.array()).max(v2.array());
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(leftBottom.array() - radius, mins);
    locateVoxelAxisIndex(rightTop.array() + radius, maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount - 1);

    edgeInds.clear();
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= surfEdgeStartInd && indI < surfTriStartInd) {
                            edgeInds.emplace_back(indI - surfEdgeStartInd);
                        }
                    }
                }
            }
        }
    }
    std::sort(edgeInds.begin(), edgeInds.end());
    edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
}

void SpatialHash::queryEdgeForEdgesWithBBoxCheck(const mesh3D& mesh,
    const vector<Eigen::Vector3d>& searchDir, double curMaxStepSize,
    int seI, std::vector<int>& sEdgeInds)
{
    if (backend == BroadPhaseBackend::LinearBVH) {
        bvh.querySweptEdgeForEdges(seI, sweptQueryRadius, sEdgeInds);
        sEdgeInds.erase(
            std::remove_if(
                sEdgeInds.begin(),
                sEdgeInds.end(),
                [seI](int candidate) { return candidate <= seI; }),
            sEdgeInds.end());
        return;
    }
    const Eigen::Matrix<double, 1, 3>& eI_v0 = mesh.vertexes[mesh.surfEdges[seI].first];
    const Eigen::Matrix<double, 1, 3>& eI_v1 = mesh.vertexes[mesh.surfEdges[seI].second];
    Eigen::Matrix<double, 1, 3> eI_v0t = eI_v0 - curMaxStepSize * searchDir[mesh.surfEdges[seI].first].transpose();
    Eigen::Matrix<double, 1, 3> eI_v1t = eI_v1 - curMaxStepSize * searchDir[mesh.surfEdges[seI].second].transpose();
    Eigen::Array<double, 1, 3> bboxEITopRight = eI_v0.array().max(eI_v0t.array()).max(eI_v1.array()).max(eI_v1t.array());
    Eigen::Array<double, 1, 3> bboxEIBottomLeft = eI_v0.array().min(eI_v0t.array()).min(eI_v1.array()).min(eI_v1t.array());
    sEdgeInds.clear();
    for (const auto& voxelInd : pointAndEdgeOccupancy[seI + surfEdgeStartInd]) {
        const auto& voxelI = voxel.find(voxelInd);
        assert(voxelI != voxel.end());
        for (const auto& indI : voxelI->second) {
            if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > seI) {
                int seJ = indI - surfEdgeStartInd;
                const Eigen::Matrix<double, 1, 3>& eJ_v0 = mesh.vertexes[mesh.surfEdges[seJ].first];
                const Eigen::Matrix<double, 1, 3>& eJ_v1 = mesh.vertexes[mesh.surfEdges[seJ].second];
                Eigen::Matrix<double, 1, 3> eJ_v0t = eJ_v0 - curMaxStepSize * searchDir[mesh.surfEdges[seJ].first].transpose();
                Eigen::Matrix<double, 1, 3> eJ_v1t = eJ_v1 - curMaxStepSize * searchDir[mesh.surfEdges[seJ].second].transpose();
                Eigen::Array<double, 1, 3> bboxEJTopRight = eJ_v0.array().max(eJ_v0t.array()).max(eJ_v1.array()).max(eJ_v1t.array());
                Eigen::Array<double, 1, 3> bboxEJBottomLeft = eJ_v0.array().min(eJ_v0t.array()).min(eJ_v1.array()).min(eJ_v1t.array());
                if (!((bboxEJBottomLeft - bboxEITopRight > 0.0).any() || (bboxEIBottomLeft - bboxEJTopRight > 0.0).any())) {
                    sEdgeInds.emplace_back(indI - surfEdgeStartInd);
                }
            }
        }
    }
    std::sort(sEdgeInds.begin(), sEdgeInds.end());
    sEdgeInds.erase(std::unique(sEdgeInds.begin(), sEdgeInds.end()), sEdgeInds.end());
}


void SpatialHash::queryEdgeForEdgesWithBBoxCheck(
    const mesh3D& mesh,
    const Eigen::Matrix<double, 1, 3>& vBegin,
    const Eigen::Matrix<double, 1, 3>& vEnd,
    double radius, std::vector<int>& edgeInds,
    int eIq)
{
    if (backend == BroadPhaseBackend::LinearBVH) {
        bvh.queryEdgeForEdges(makeEdgeAABB(vBegin.transpose(), vEnd.transpose()), radius, edgeInds);
        edgeInds.erase(
            std::remove_if(
                edgeInds.begin(),
                edgeInds.end(),
                [eIq](int candidate) { return candidate <= eIq; }),
            edgeInds.end());
        return;
    }
    // timer_mt.start(19);
    Eigen::Matrix<double, 1, 3> leftBottom = vBegin.array().min(vEnd.array()) - radius;
    Eigen::Matrix<double, 1, 3> rightTop = vBegin.array().max(vEnd.array()) + radius;
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(leftBottom, mins);
    locateVoxelAxisIndex(rightTop, maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount.array() - 1);

    // timer_mt.start(20);
    edgeInds.resize(0);
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                // timer_mt.start(21);
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > eIq) {
                            int seJ = indI - surfEdgeStartInd;
                            const Eigen::Matrix<double, 1, 3>& eJ_v0 = mesh.vertexes[mesh.surfEdges[seJ].first];
                            const Eigen::Matrix<double, 1, 3>& eJ_v1 = mesh.vertexes[mesh.surfEdges[seJ].second];
                            Eigen::Array<double, 1, 3> bboxEJTopRight = eJ_v0.array().max(eJ_v1.array());
                            Eigen::Array<double, 1, 3> bboxEJBottomLeft = eJ_v0.array().min(eJ_v1.array());
                            if (!((bboxEJBottomLeft - rightTop.array() > 0.0).any() || (leftBottom.array() - bboxEJTopRight > 0.0).any())) {
                                edgeInds.emplace_back(seJ);
                            }
                        }
                    }
                }
                // timer_mt.start(20);
            }
        }
    }
    std::sort(edgeInds.begin(), edgeInds.end());
    edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
    // timer_mt.stop();
}

void SpatialHash::queryEdgeForPE(const Vector3d& vBegin, const Vector3d& vEnd, std::vector<int>& svInds, std::vector<int>& edgeInds)
{
    Vector3d leftBottom = vBegin.array().min(vEnd.array());
    Vector3d rightTop = vBegin.array().max(vEnd.array());
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(leftBottom, mins);
    locateVoxelAxisIndex(rightTop, maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount - 1);

    // timer_mt.start(20);
    svInds.resize(0);
    edgeInds.resize(0);
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                // timer_mt.start(21);
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI < surfEdgeStartInd) {
                            svInds.emplace_back(indI);
                        }
                        else if (indI < surfTriStartInd) {
                            edgeInds.emplace_back(indI - surfEdgeStartInd);
                        }
                    }
                }
                // timer_mt.start(20);
            }
        }
    }
    std::sort(edgeInds.begin(), edgeInds.end());
    edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
    std::sort(svInds.begin(), svInds.end());
    svInds.erase(std::unique(svInds.begin(), svInds.end()), svInds.end());
}

void SpatialHash::queryPointForPrimitives(const Vector3d& pos, const Vector3d& dir, std::unordered_set<int>& sVInds, std::unordered_set<int>& sEdgeInds, std::unordered_set<int>& sTriInds)
{
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(pos.array().min((pos + dir).array()), mins);
    locateVoxelAxisIndex(pos.array().max((pos + dir).array()), maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount - 1);

    sVInds.clear();
    sEdgeInds.clear();
    sTriInds.clear();
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI < surfEdgeStartInd) {
                            sVInds.insert(indI);
                        }
                        else if (indI < surfTriStartInd) {
                            sEdgeInds.insert(indI - surfEdgeStartInd);
                        }
                        else {
                            sTriInds.insert(indI - surfTriStartInd);
                        }
                    }
                }
            }
        }
    }
}

void SpatialHash::calculateActivateSet(mesh3D& mesh) {
    //mesh.Self_ActiveSet.resize(0);
    double sqrtDHat = std::sqrt(mesh.Hhat);
    vector<vector<EncodedContact>> constraintSetPT(mesh.surfVerts.size());
    vector<vector<int>> partialCCDSetPT(mesh.surfVerts.size());
    tbb::enumerable_thread_specific<std::vector<int>> ptCandidateStorage;

    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI)
        //for (int svI = 0; svI < mesh.surfVerts.size(); ++svI)
        {
            int vI = mesh.surfVerts[svI];
            std::vector<int>& triInds = ptCandidateStorage.local();
            queryPointForTriangles(mesh.vertexes[vI], sqrtDHat, triInds);
            const CPUAABB pointBox = makePointAABB(mesh.vertexes[vI]);
            for (const auto& sfI : triInds)
            {
                const Vector4i& sfVInd = mesh.surface[sfI];
                if (isExternalColliderBoundary(mesh.boundaryTypes[vI])
                    && isExternalColliderBoundary(mesh.boundaryTypes[sfVInd[0]])
                    && isExternalColliderBoundary(mesh.boundaryTypes[sfVInd[1]])
                    && isExternalColliderBoundary(mesh.boundaryTypes[sfVInd[2]])) {
                    continue;
                }
                if (vI == sfVInd[0] || vI == sfVInd[1] || vI == sfVInd[2]) {
                    continue;
                }
                const CPUAABB triangleBox = makeTriangleAABB(
                    mesh.vertexes[sfVInd[0]],
                    mesh.vertexes[sfVInd[1]],
                    mesh.vertexes[sfVInd[2]]);
                if (!pointBox.overlaps(triangleBox, sqrtDHat)) {
                    continue;
                }
                partialCCDSetPT[svI].emplace_back(sfI);

                    int dtype = dType_PT(mesh.vertexes[vI], mesh.vertexes[sfVInd[0]], mesh.vertexes[sfVInd[1]], mesh.vertexes[sfVInd[2]]);
                    double d;
                    switch (dtype) {
                    case 0: {
                        d_PP(mesh.vertexes[vI], mesh.vertexes[sfVInd[0]], d);
                        if (d < mesh.Hhat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[0], -1, -1);
                        }
                        break;
                    }

                    case 1: {
                        d_PP(mesh.vertexes[vI], mesh.vertexes[sfVInd[1]], d);
                        if (d < mesh.Hhat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[1], -1, -1);
                        }
                        break;
                    }

                    case 2: {
                        d_PP(mesh.vertexes[vI], mesh.vertexes[sfVInd[2]], d);
                        if (d < mesh.Hhat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[2], -1, -1);
                        }
                        break;
                    }

                    case 3: {
                        d_PE(mesh.vertexes[vI], mesh.vertexes[sfVInd[0]], mesh.vertexes[sfVInd[1]], d);
                        if (d < mesh.Hhat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[0], sfVInd[1], -1);
                        }
                        break;
                    }

                    case 4: {
                        d_PE(mesh.vertexes[vI], mesh.vertexes[sfVInd[1]], mesh.vertexes[sfVInd[2]], d);
                        if (d < mesh.Hhat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[1], sfVInd[2], -1);
                        }
                        break;
                    }

                    case 5: {
                        d_PE(mesh.vertexes[vI], mesh.vertexes[sfVInd[2]], mesh.vertexes[sfVInd[0]], d);
                        if (d < mesh.Hhat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[2], sfVInd[0], -1);
                        }
                        break;
                    }

                    case 6: {
                        d_PT(mesh.vertexes[vI], mesh.vertexes[sfVInd[0]], mesh.vertexes[sfVInd[1]], mesh.vertexes[sfVInd[2]], d);
                        if (d < mesh.Hhat) {
                            constraintSetPT[svI].emplace_back(-vI - 1, sfVInd[0], sfVInd[1], sfVInd[2]);
                        }
                        break;
                    }

                    default:
                        break;
                    }
            }
        }
    );
    vector<vector<EncodedContact>> constraintSetEE(mesh.surfEdges.size());
    vector<vector<int>> partialCCDSetEE(mesh.surfEdges.size());
    tbb::enumerable_thread_specific<std::vector<int>> eeCandidateStorage;

    tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int eI)
    //for (int eI = 0; eI < mesh.surfEdges.size(); ++eI)
    {
        const auto& meshEI = mesh.surfEdges[eI];
        vector<int>& edgeInds = eeCandidateStorage.local();
        queryEdgeForEdgesWithBBoxCheck(mesh, mesh.vertexes[meshEI.first].transpose(), mesh.vertexes[meshEI.second].transpose(), sqrtDHat, edgeInds, eI);
        const CPUAABB firstEdgeBox = makeEdgeAABB(
            mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second]);
        for (const auto& eJ : edgeInds) {
            const auto& meshEJ = mesh.surfEdges[eJ];
            if (isExternalColliderBoundary(mesh.boundaryTypes[meshEI.first])
                && isExternalColliderBoundary(mesh.boundaryTypes[meshEI.second])
                && isExternalColliderBoundary(mesh.boundaryTypes[meshEJ.first])
                && isExternalColliderBoundary(mesh.boundaryTypes[meshEJ.second])) {
                continue;
            }
            if (meshEI.first == meshEJ.first || meshEI.first == meshEJ.second
                || meshEI.second == meshEJ.first || meshEI.second == meshEJ.second
                || eI > eJ) {
                continue;
            }
            const CPUAABB secondEdgeBox = makeEdgeAABB(
                mesh.vertexes[meshEJ.first], mesh.vertexes[meshEJ.second]);
            if (!firstEdgeBox.overlaps(secondEdgeBox, sqrtDHat)) {
                continue;
            }
            partialCCDSetEE[eI].emplace_back(eJ);

                int dtype = dType_EE(mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second], mesh.vertexes[meshEJ.first], mesh.vertexes[meshEJ.second]);

                double EECrossSqNorm, eps_x;
                computeEECrossSqNorm(mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second], mesh.vertexes[meshEJ.first], mesh.vertexes[meshEJ.second], EECrossSqNorm);
                compute_eps_x(mesh, meshEI.first, meshEI.second, meshEJ.first, meshEJ.second, eps_x);

                int add_e = (EECrossSqNorm < eps_x) ? -eJ - 2 : -1;
                double d;
                switch (dtype) {
                case 0: {
                    d_PP(mesh.vertexes[meshEI.first], mesh.vertexes[meshEJ.first], d);
                    if (d < mesh.Hhat) {
                        constraintSetEE[eI].emplace_back(-meshEI.first - 1, meshEJ.first, -1, add_e);
                    }
                    break;
                }

                case 1: {
                    d_PP(mesh.vertexes[meshEI.first], mesh.vertexes[meshEJ.second], d);
                    if (d < mesh.Hhat) {
                        constraintSetEE[eI].emplace_back(-meshEI.first - 1, meshEJ.second, -1, add_e);
                    }
                    break;
                }

                case 2: {
                    d_PE(mesh.vertexes[meshEI.first], mesh.vertexes[meshEJ.first], mesh.vertexes[meshEJ.second], d);
                    if (d < mesh.Hhat) {
                        constraintSetEE[eI].emplace_back(-meshEI.first - 1, meshEJ.first, meshEJ.second, add_e);
                    }
                    break;
                }

                case 3: {
                    d_PP(mesh.vertexes[meshEI.second], mesh.vertexes[meshEJ.first], d);
                    if (d < mesh.Hhat) {
                        constraintSetEE[eI].emplace_back(-meshEI.second - 1, meshEJ.first, -1, add_e);
                    }
                    break;
                }

                case 4: {
                    d_PP(mesh.vertexes[meshEI.second], mesh.vertexes[meshEJ.second], d);
                    if (d < mesh.Hhat) {
                        constraintSetEE[eI].emplace_back(-meshEI.second - 1, meshEJ.second, -1, add_e);
                    }
                    break;
                }

                case 5: {
                    d_PE(mesh.vertexes[meshEI.second], mesh.vertexes[meshEJ.first], mesh.vertexes[meshEJ.second], d);
                    if (d < mesh.Hhat) {
                        constraintSetEE[eI].emplace_back(-meshEI.second - 1, meshEJ.first, meshEJ.second, add_e);
                    }
                    break;
                }

                case 6: {
                    d_PE(mesh.vertexes[meshEJ.first], mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second], d);
                    if (d < mesh.Hhat) {
                        constraintSetEE[eI].emplace_back(-meshEJ.first - 1, meshEI.first, meshEI.second, add_e);
                    }
                    break;
                }

                case 7: {
                    d_PE(mesh.vertexes[meshEJ.second], mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second], d);
                    if (d < mesh.Hhat) {
                        constraintSetEE[eI].emplace_back(-meshEJ.second - 1, meshEI.first, meshEI.second, add_e);
                    }
                    break;
                }

                case 8: {
                    d_EE(mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second], mesh.vertexes[meshEJ.first], mesh.vertexes[meshEJ.second], d);
                    if (d < mesh.Hhat) {
                        if (add_e <= -2) {
                            constraintSetEE[eI].emplace_back(meshEI.first, meshEI.second, meshEJ.first, -meshEJ.second - mesh.surfEdges.size() - 2);
                        }
                        else {
                            constraintSetEE[eI].emplace_back(meshEI.first, meshEI.second, meshEJ.first, meshEJ.second);
                        }
                    }
                    break;
                }

                default:
                    break;
                }
        }
    }
    );

    size_t partialCCDPairCount = 0;
    for (const auto& pairs : partialCCDSetPT) {
        partialCCDPairCount += pairs.size();
    }
    for (const auto& pairs : partialCCDSetEE) {
        partialCCDPairCount += pairs.size();
    }
    mesh.Self_CCD_ActiveSet.clear();
    mesh.Self_CCD_ActiveSet.reserve(partialCCDPairCount);
    for (int surfaceVertex = 0; surfaceVertex < static_cast<int>(partialCCDSetPT.size()); ++surfaceVertex) {
        for (const int face : partialCCDSetPT[surfaceVertex]) {
            mesh.Self_CCD_ActiveSet.emplace_back(-surfaceVertex - 1, face);
        }
    }
    for (int firstEdge = 0; firstEdge < static_cast<int>(partialCCDSetEE.size()); ++firstEdge) {
        for (const int secondEdge : partialCCDSetEE[firstEdge]) {
            mesh.Self_CCD_ActiveSet.emplace_back(firstEdge, secondEdge);
        }
    }

    mesh.Self_ActiveSet.resize(0);
    mesh.Self_ActiveSet.reserve(constraintSetPT.size() + constraintSetEE.size());

    std::map<EncodedContact, int> constraintCounter;
    for (const auto& csI : constraintSetPT) {
        for (const auto& cI : csI) {
            if (cI[3] < 0) {
                // PP or PE
                ++constraintCounter[cI];
            }
            else {
                mesh.Self_ActiveSet.emplace_back(cI);
            }
        }
    }
    mesh.Self_EE_ActiveSet.resize(0);
    mesh.Self_EEeIe_ActiveSet.resize(0);
    int eI = 0;
    for (const auto& csI : constraintSetEE) {
        for (const auto& cI : csI) {
            if (cI[3] >= 0) {
                // regular EE
                mesh.Self_ActiveSet.emplace_back(cI);
            }
            else if (cI[3] == -1) {
                // regular PP or PE
                ++constraintCounter[cI];
            }
            else if (cI[3] >= -mesh.surfEdges.size() - 1) {
                // nearly parallel PP or PE
                mesh.Self_EE_ActiveSet.emplace_back(cI[0], cI[1], cI[2], -1);
                mesh.Self_EEeIe_ActiveSet.emplace_back(eI, -cI[3] - 2);
            }
            else {
                // nearly parallel EE
                mesh.Self_EE_ActiveSet.emplace_back(
                    cI[0], cI[1], cI[2],
                    -cI[3] - static_cast<int>(mesh.surfEdges.size()) - 2);
                mesh.Self_EEeIe_ActiveSet.emplace_back(-1, -1);
            }
        }
        ++eI;
    }

    mesh.Self_ActiveSet.reserve(mesh.Self_ActiveSet.size() + constraintCounter.size());
    for (const auto& ccI : constraintCounter) {
        mesh.Self_ActiveSet.emplace_back(
            ccI.first[0], ccI.first[1], ccI.first[2], -ccI.second);
    }
}


void Ground::init(const Vector3d& m_normal, const double& d) {
    normal = m_normal;
    D = d;
}

Ground::Ground() {
    init(Vector3d(0, 1, 0), -1);
}

double Ground::calculateGapFromObj(const mesh3D& mesh, const int& vId) const {
    double dist = normal.dot(mesh.vertexes[vId]) - D;
    return dist * dist;
}

void Ground::calculateActivateSet(mesh3D& mesh) const {
    mesh.Environment_ActiveSet.resize(0);


    tbb::spin_mutex groundMutex;//, countMutex3;
    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int i)
        //for (int i = 0;i < mesh.surfVerts.size();i++) 
        {
            double dis = calculateGapFromObj(mesh, mesh.surfVerts[i]);
            if (dis < mesh.Hhat) {
                groundMutex.lock();
                mesh.Environment_ActiveSet.push_back(mesh.surfVerts[i]);
                groundMutex.unlock();
            }
        }
    );
}
