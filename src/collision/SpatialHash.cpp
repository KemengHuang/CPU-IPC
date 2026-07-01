#include "math/Constants.h"
#include "core/IPCDistance.h"
#include "core/SimulationParameters.h"
#include "collision/SpatialHash.h"
#include <algorithm>
#include <assert.h>
#include <map>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>

using namespace std;
using namespace Eigen;

void SpatialHash::locateVoxelAxisIndex(const Vector3d& pos, Eigen::Array<int, 1, 3>& voxelAxisIndex) const
{
    voxelAxisIndex = ((pos - leftBottomCorner) * one_div_voxelSize).array().floor().template cast<int>();
}

int SpatialHash::locateVoxelIndex(const Vector3d& pos) const
{
    Eigen::Array<int, 1, 3> voxelAxisIndex;
    locateVoxelAxisIndex(pos, voxelAxisIndex);
    return voxelAxisIndex2VoxelIndex(voxelAxisIndex.data());
}

int SpatialHash::voxelAxisIndex2VoxelIndex(const int voxelAxisIndex[3]) const
{
    return voxelAxisIndex2VoxelIndex(voxelAxisIndex[0], voxelAxisIndex[1], voxelAxisIndex[2]);
}

int SpatialHash::voxelAxisIndex2VoxelIndex(int ix, int iy, int iz) const
{
    return ix + iy * voxelCount[0] + iz * voxelCount0x1;
}

void SpatialHash::computeVoxelGrid(const mesh3D& mesh, double voxelSize)
{
    one_div_voxelSize = 1.0 / voxelSize;
    Eigen::Array<double, 1, 3> range = rightTopCorner - leftBottomCorner;
    voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
    if (voxelCount.minCoeff() <= 0) {
        // cast overflow due to huge search direction
        one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
        voxelCount.setOnes();
    }
    voxelCount0x1 = voxelCount[0] * voxelCount[1];
}

void SpatialHash::build(const mesh3D& mesh, double voxelSize)
{
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
        });

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
        });

    computeVoxelGrid(mesh, voxelSize);

    surfEdgeStartInd = mesh.surfVerts.size();
    surfTriStartInd = surfEdgeStartInd + mesh.surfEdges.size();

    std::vector<Eigen::Array<int, 1, 3>> svVoxelAxisIndex(mesh.surfVerts.size());
    std::vector<int> vI2SVI(mesh.vertexNum);

    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI) {
        int vI = mesh.surfVerts[svI];
        locateVoxelAxisIndex(mesh.vertexes[vI], svVoxelAxisIndex[svI]);
        vI2SVI[vI] = svI;
    });

    voxel.clear();

    for (int svI = 0; svI < mesh.surfVerts.size(); ++svI) {
        voxel[locateVoxelIndex(mesh.vertexes[mesh.surfVerts[svI]])].emplace_back(svI);
    }

    std::vector<std::vector<int>> voxelLoc_e(mesh.surfEdges.size());

    tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int seCount) {
        const auto& seI = mesh.surfEdges[seCount];

        const Eigen::Array<int, 1, 3>& voxelAxisIndex_first = svVoxelAxisIndex[vI2SVI[seI.first]];
        const Eigen::Array<int, 1, 3>& voxelAxisIndex_second = svVoxelAxisIndex[vI2SVI[seI.second]];
        Eigen::Array<int, 1, 3> mins = voxelAxisIndex_first.min(voxelAxisIndex_second);
        Eigen::Array<int, 1, 3> maxs = voxelAxisIndex_first.max(voxelAxisIndex_second);
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    voxelLoc_e[seCount].emplace_back(ix + yzOffset);
                }
            }
        }
    });

    std::vector<std::vector<int>> voxelLoc_sf(mesh.surface.size());

    tbb::parallel_for(0, (int)mesh.surface.size(), 1, [&](int sfI) {
        const Eigen::Array<int, 1, 3>& voxelAxisIndex0 = svVoxelAxisIndex[vI2SVI[mesh.surface[sfI][0]]];
        const Eigen::Array<int, 1, 3>& voxelAxisIndex1 = svVoxelAxisIndex[vI2SVI[mesh.surface[sfI][1]]];
        const Eigen::Array<int, 1, 3>& voxelAxisIndex2 = svVoxelAxisIndex[vI2SVI[mesh.surface[sfI][2]]];
        Eigen::Array<int, 1, 3> mins = voxelAxisIndex0.min(voxelAxisIndex1).min(voxelAxisIndex2);
        Eigen::Array<int, 1, 3> maxs = voxelAxisIndex0.max(voxelAxisIndex1).max(voxelAxisIndex2);
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    voxelLoc_sf[sfI].emplace_back(ix + yzOffset);
                }
            }
        }
    });

    for (int seCount = 0; seCount < voxelLoc_e.size(); ++seCount) {
        for (const auto& voxelI : voxelLoc_e[seCount]) {
            voxel[voxelI].emplace_back(seCount + surfEdgeStartInd);
        }
    }
    for (int sfI = 0; sfI < voxelLoc_sf.size(); ++sfI) {
        for (const auto& voxelI : voxelLoc_sf[sfI]) {
            voxel[voxelI].emplace_back(sfI + surfTriStartInd);
        }
    }
}

void SpatialHash::build(const mesh3D& mesh, const vector<Vector3d>& searchDir, double& curMaxStepSize, double voxelSize, bool use_V_prev)
{
    double pSize = 0;
    for (int svI = 0; svI < mesh.surfVerts.size(); ++svI) {
        int vI = mesh.surfVerts[svI];
        pSize += std::abs(searchDir[vI][0]);
        pSize += std::abs(searchDir[vI][1]);
        pSize += std::abs(searchDir[vI][2]);
    }
    pSize /= mesh.surfVerts.size() * 3;

    const double spanSize = curMaxStepSize * pSize / voxelSize;

    const vector<Vector3d>& V = use_V_prev ? mesh.V_prev : mesh.vertexes;
    vector<Vector3d> SVt(mesh.surfVerts.size());
    std::vector<int> vI2SVI(mesh.vertexNum);

    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI) {
        int vI = mesh.surfVerts[svI];
        SVt[svI] = V[vI] - curMaxStepSize * searchDir[vI];
        vI2SVI[vI] = svI;
    });

    leftBottomCorner = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, (int)mesh.surfVerts.size()), Vector3d(1e32, 1e32, 1e32), [&](const tbb::blocked_range<int>& rg, Vector3d sceneSize) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                Vector3d pos = V[mesh.surfVerts[i]];
                sceneSize[0] = min(pos[0], sceneSize[0]);
                sceneSize[1] = min(pos[1], sceneSize[1]);
                sceneSize[2] = min(pos[2], sceneSize[2]);

                pos = SVt[i];
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
        });

    rightTopCorner = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, (int)mesh.surfVerts.size()), Vector3d(-1e32, -1e32, -1e32), [&](const tbb::blocked_range<int>& rg, Vector3d sceneSize) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                Vector3d pos = V[mesh.surfVerts[i]];
                sceneSize[0] = max(pos[0], sceneSize[0]);
                sceneSize[1] = max(pos[1], sceneSize[1]);
                sceneSize[2] = max(pos[2], sceneSize[2]);

                pos = SVt[i];
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
        });

    computeVoxelGrid(mesh, voxelSize);

    surfEdgeStartInd = mesh.surfVerts.size();
    surfTriStartInd = surfEdgeStartInd + mesh.surfEdges.size();

    std::vector<Eigen::Array<int, 1, 3>> svMinVAI(mesh.surfVerts.size());
    std::vector<Eigen::Array<int, 1, 3>> svMaxVAI(mesh.surfVerts.size());

    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI) {
        int vI = mesh.surfVerts[svI];
        Eigen::Array<int, 1, 3> v0VAI, vtVAI;
        locateVoxelAxisIndex(V[vI], v0VAI);
        locateVoxelAxisIndex(SVt[svI], vtVAI);
        svMinVAI[svI] = v0VAI.min(vtVAI);
        svMaxVAI[svI] = v0VAI.max(vtVAI);
    });

    voxel.clear();
    pointAndEdgeOccupancy.resize(0);
    pointAndEdgeOccupancy.resize(surfTriStartInd);

    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI) {
        const Eigen::Array<int, 1, 3>& mins = svMinVAI[svI];
        const Eigen::Array<int, 1, 3>& maxs = svMaxVAI[svI];
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
    });

    tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int seCount) {
        int seIInd = seCount + surfEdgeStartInd;
        const auto& seI = mesh.surfEdges[seCount];

        Eigen::Array<int, 1, 3> mins = svMinVAI[vI2SVI[seI.first]].min(svMinVAI[vI2SVI[seI.second]]);
        Eigen::Array<int, 1, 3> maxs = svMaxVAI[vI2SVI[seI.first]].max(svMaxVAI[vI2SVI[seI.second]]);
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
    });

    std::vector<std::vector<int>> voxelLoc_sf(mesh.surface.size());

    tbb::parallel_for(0, (int)mesh.surface.size(), 1, [&](int sfI) {
        Eigen::Array<int, 1, 3> mins = svMinVAI[vI2SVI[mesh.surface[sfI][0]]].min(svMinVAI[vI2SVI[mesh.surface[sfI][1]]]).min(svMinVAI[vI2SVI[mesh.surface[sfI][2]]]);
        Eigen::Array<int, 1, 3> maxs = svMaxVAI[vI2SVI[mesh.surface[sfI][0]]].max(svMaxVAI[vI2SVI[mesh.surface[sfI][1]]]).max(svMaxVAI[vI2SVI[mesh.surface[sfI][2]]]);
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    voxelLoc_sf[sfI].emplace_back(ix + yzOffset);
                }
            }
        }
    });

    for (int i = 0; i < pointAndEdgeOccupancy.size(); ++i) {
        for (const auto& voxelI : pointAndEdgeOccupancy[i]) {
            voxel[voxelI].emplace_back(i);
        }
    }
    for (int sfI = 0; sfI < voxelLoc_sf.size(); ++sfI) {
        for (const auto& voxelI : voxelLoc_sf[sfI]) {
            voxel[voxelI].emplace_back(sfI + surfTriStartInd);
        }
    }
}


void SpatialHash::queryPointForTriangles(const Vector3d& pos, double radius, std::unordered_set<int>& triInds) const
{
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(pos.array() - radius, mins);
    locateVoxelAxisIndex(pos.array() + radius, maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount.array() - 1);

    thread_local std::vector<int> buffer;
    buffer.clear();
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= surfTriStartInd) {
                            buffer.emplace_back(indI - surfTriStartInd);
                        }
                    }
                }
            }
        }
    }
    std::sort(buffer.begin(), buffer.end());
    buffer.erase(std::unique(buffer.begin(), buffer.end()), buffer.end());

    triInds.clear();
    triInds.insert(buffer.begin(), buffer.end());
}

void SpatialHash::queryPointForTriangles(int svI, std::unordered_set<int>& sTriInds) const
{
    thread_local std::vector<int> buffer;
    buffer.clear();
    for (const auto& voxelInd : pointAndEdgeOccupancy[svI]) {
        const auto& voxelI = voxel.find(voxelInd);
        assert(voxelI != voxel.end());
        for (const auto& indI : voxelI->second) {
            if (indI >= surfTriStartInd) {
                buffer.emplace_back(indI - surfTriStartInd);
            }
        }
    }
    std::sort(buffer.begin(), buffer.end());
    buffer.erase(std::unique(buffer.begin(), buffer.end()), buffer.end());

    sTriInds.clear();
    sTriInds.insert(buffer.begin(), buffer.end());
}

void SpatialHash::queryTriangleForPoints(const Vector3d& v0,
    const Vector3d& v1,
    const Vector3d& v2,
    double radius, std::unordered_set<int>& pointInds) const
{
    Vector3d leftBottom = v0.array().min(v1.array()).min(v2.array());
    Vector3d rightTop = v0.array().max(v1.array()).max(v2.array());
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(leftBottom.array() - radius, mins);
    locateVoxelAxisIndex(rightTop.array() + radius, maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount - 1);

    thread_local std::vector<int> buffer;
    buffer.clear();
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI < surfEdgeStartInd) {
                            buffer.emplace_back(indI);
                        }
                    }
                }
            }
        }
    }
    std::sort(buffer.begin(), buffer.end());
    buffer.erase(std::unique(buffer.begin(), buffer.end()), buffer.end());

    pointInds.clear();
    pointInds.insert(buffer.begin(), buffer.end());
}

void SpatialHash::queryTriangleForEdges(const Vector3d& v0,
    const Vector3d& v1,
    const Vector3d& v2,
    double radius, std::unordered_set<int>& edgeInds) const
{
    Vector3d leftBottom = v0.array().min(v1.array()).min(v2.array());
    Vector3d rightTop = v0.array().max(v1.array()).max(v2.array());
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(leftBottom.array() - radius, mins);
    locateVoxelAxisIndex(rightTop.array() + radius, maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount - 1);

    thread_local std::vector<int> buffer;
    buffer.clear();
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI >= surfEdgeStartInd && indI < surfTriStartInd) {
                            buffer.emplace_back(indI - surfEdgeStartInd);
                        }
                    }
                }
            }
        }
    }
    std::sort(buffer.begin(), buffer.end());
    buffer.erase(std::unique(buffer.begin(), buffer.end()), buffer.end());

    edgeInds.clear();
    edgeInds.insert(buffer.begin(), buffer.end());
}

void SpatialHash::queryEdgeForEdgesWithBBoxCheck(const mesh3D& mesh,
    const vector<Eigen::Vector3d>& searchDir, double curMaxStepSize,
    int seI, std::unordered_set<int>& sEdgeInds) const
{
    const Eigen::Matrix<double, 1, 3>& eI_v0 = mesh.vertexes[mesh.surfEdges[seI].first];
    const Eigen::Matrix<double, 1, 3>& eI_v1 = mesh.vertexes[mesh.surfEdges[seI].second];
    Eigen::Matrix<double, 1, 3> eI_v0t = eI_v0 - curMaxStepSize * searchDir[mesh.surfEdges[seI].first].transpose();
    Eigen::Matrix<double, 1, 3> eI_v1t = eI_v1 - curMaxStepSize * searchDir[mesh.surfEdges[seI].second].transpose();
    Eigen::Array<double, 1, 3> bboxEITopRight = eI_v0.array().max(eI_v0t.array()).max(eI_v1.array()).max(eI_v1t.array());
    Eigen::Array<double, 1, 3> bboxEIBottomLeft = eI_v0.array().min(eI_v0t.array()).min(eI_v1.array()).min(eI_v1t.array());

    thread_local std::vector<int> buffer;
    buffer.clear();
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
                    buffer.emplace_back(indI - surfEdgeStartInd);
                }
            }
        }
    }
    std::sort(buffer.begin(), buffer.end());
    buffer.erase(std::unique(buffer.begin(), buffer.end()), buffer.end());

    sEdgeInds.clear();
    sEdgeInds.insert(buffer.begin(), buffer.end());
}

void SpatialHash::queryEdgeForEdgesWithBBoxCheck(
    const mesh3D& mesh,
    const Eigen::Matrix<double, 1, 3>& vBegin,
    const Eigen::Matrix<double, 1, 3>& vEnd,
    double radius, std::vector<int>& edgeInds,
    int eIq) const
{
    Eigen::Matrix<double, 1, 3> leftBottom = vBegin.array().min(vEnd.array()) - radius;
    Eigen::Matrix<double, 1, 3> rightTop = vBegin.array().max(vEnd.array()) + radius;
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(leftBottom, mins);
    locateVoxelAxisIndex(rightTop, maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount.array() - 1);

    edgeInds.resize(0);
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
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
            }
        }
    }
    std::sort(edgeInds.begin(), edgeInds.end());
    edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
}

void SpatialHash::queryEdgeForPE(const Vector3d& vBegin, const Vector3d& vEnd, std::vector<int>& svInds, std::vector<int>& edgeInds) const
{
    Vector3d leftBottom = vBegin.array().min(vEnd.array());
    Vector3d rightTop = vBegin.array().max(vEnd.array());
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(leftBottom, mins);
    locateVoxelAxisIndex(rightTop, maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount - 1);

    svInds.resize(0);
    edgeInds.resize(0);
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
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
            }
        }
    }
    std::sort(edgeInds.begin(), edgeInds.end());
    edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
    std::sort(svInds.begin(), svInds.end());
    svInds.erase(std::unique(svInds.begin(), svInds.end()), svInds.end());
}

void SpatialHash::queryPointForPrimitives(const Vector3d& pos, const Vector3d& dir, std::unordered_set<int>& sVInds, std::unordered_set<int>& sEdgeInds, std::unordered_set<int>& sTriInds) const
{
    Eigen::Array<int, 1, 3> mins, maxs;
    locateVoxelAxisIndex(pos.array().min((pos + dir).array()), mins);
    locateVoxelAxisIndex(pos.array().max((pos + dir).array()), maxs);
    mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
    maxs = maxs.min(voxelCount - 1);

    thread_local std::vector<int> vBuffer, eBuffer, tBuffer;
    vBuffer.clear();
    eBuffer.clear();
    tBuffer.clear();
    for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
        int zOffset = iz * voxelCount0x1;
        for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
            int yzOffset = iy * voxelCount[0] + zOffset;
            for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                const auto voxelI = voxel.find(ix + yzOffset);
                if (voxelI != voxel.end()) {
                    for (const auto& indI : voxelI->second) {
                        if (indI < surfEdgeStartInd) {
                            vBuffer.emplace_back(indI);
                        }
                        else if (indI < surfTriStartInd) {
                            eBuffer.emplace_back(indI - surfEdgeStartInd);
                        }
                        else {
                            tBuffer.emplace_back(indI - surfTriStartInd);
                        }
                    }
                }
            }
        }
    }
    std::sort(vBuffer.begin(), vBuffer.end());
    vBuffer.erase(std::unique(vBuffer.begin(), vBuffer.end()), vBuffer.end());
    std::sort(eBuffer.begin(), eBuffer.end());
    eBuffer.erase(std::unique(eBuffer.begin(), eBuffer.end()), eBuffer.end());
    std::sort(tBuffer.begin(), tBuffer.end());
    tBuffer.erase(std::unique(tBuffer.begin(), tBuffer.end()), tBuffer.end());

    sVInds.clear();
    sVInds.insert(vBuffer.begin(), vBuffer.end());
    sEdgeInds.clear();
    sEdgeInds.insert(eBuffer.begin(), eBuffer.end());
    sTriInds.clear();
    sTriInds.insert(tBuffer.begin(), tBuffer.end());
}


void SpatialHash::collectPointTrianglePairs(mesh3D& mesh, std::vector<MMCVID>& ptPairs)
{
    const double sqrtDHat = std::sqrt(mesh.Hhat);
    std::vector<std::vector<MMCVID>> constraintSetPT(mesh.surfVerts.size());
    std::vector<std::vector<std::pair<int, int>>> ccdPairs(mesh.surfVerts.size());

    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI) {
        const int vI = mesh.surfVerts[svI];
        thread_local std::vector<int> triInds;
        triInds.clear();

        Eigen::Array<int, 1, 3> mins, maxs;
        locateVoxelAxisIndex(mesh.vertexes[vI].array() - sqrtDHat, mins);
        locateVoxelAxisIndex(mesh.vertexes[vI].array() + sqrtDHat, maxs);
        mins = mins.max(Eigen::Array<int, 1, 3>::Zero());
        maxs = maxs.min(voxelCount.array() - 1);

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

        for (const int sfI : triInds) {
            const RowVector4i sfVInd = mesh.surface[sfI].transpose();
            if (mesh.boundaryTypes[vI] >= 2 && mesh.boundaryTypes[sfVInd[0]] >= 2 && mesh.boundaryTypes[sfVInd[1]] >= 2 && mesh.boundaryTypes[sfVInd[2]] >= 2) {
                continue;
            }
            if (vI == sfVInd[0] || vI == sfVInd[1] || vI == sfVInd[2]) {
                continue;
            }

            ccdPairs[svI].emplace_back(-svI - 1, sfI);

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
    });

    size_t totalPairs = 0;
    size_t totalCCD = 0;
    for (size_t svI = 0; svI < mesh.surfVerts.size(); ++svI) {
        totalPairs += constraintSetPT[svI].size();
        totalCCD += ccdPairs[svI].size();
    }

    ptPairs.clear();
    ptPairs.reserve(totalPairs);
    mesh.Self_CCD_ActiveSet.reserve(mesh.Self_CCD_ActiveSet.size() + totalCCD);
    for (size_t svI = 0; svI < mesh.surfVerts.size(); ++svI) {
        ptPairs.insert(ptPairs.end(), constraintSetPT[svI].begin(), constraintSetPT[svI].end());
        mesh.Self_CCD_ActiveSet.insert(mesh.Self_CCD_ActiveSet.end(), ccdPairs[svI].begin(), ccdPairs[svI].end());
    }
}

void SpatialHash::collectEdgeEdgePairs(mesh3D& mesh, std::vector<std::vector<MMCVID>>& eePairsPerEdge)
{
    const double sqrtDHat = std::sqrt(mesh.Hhat);
    eePairsPerEdge.assign(mesh.surfEdges.size(), std::vector<MMCVID>());
    std::vector<std::vector<std::pair<int, int>>> ccdPairs(mesh.surfEdges.size());

    tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int eI) {
        const auto& meshEI = mesh.surfEdges[eI];
        thread_local std::vector<int> edgeInds;
        edgeInds.clear();

        queryEdgeForEdgesWithBBoxCheck(mesh, mesh.vertexes[meshEI.first].transpose(), mesh.vertexes[meshEI.second].transpose(), sqrtDHat, edgeInds, eI);

        for (const int eJ : edgeInds) {
            const auto& meshEJ = mesh.surfEdges[eJ];
            if (mesh.boundaryTypes[meshEI.first] >= 2 && mesh.boundaryTypes[meshEI.second] >= 2 && mesh.boundaryTypes[meshEJ.first] >= 2 && mesh.boundaryTypes[meshEJ.second] >= 2) {
                continue;
            }
            if (meshEI.first == meshEJ.first || meshEI.first == meshEJ.second || meshEI.second == meshEJ.first || meshEI.second == meshEJ.second || eI > eJ) {
                continue;
            }

            ccdPairs[eI].emplace_back(eI, eJ);

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
                    eePairsPerEdge[eI].emplace_back(-meshEI.first - 1, meshEJ.first, -1, add_e);
                }
                break;
            }
            case 1: {
                d_PP(mesh.vertexes[meshEI.first], mesh.vertexes[meshEJ.second], d);
                if (d < mesh.Hhat) {
                    eePairsPerEdge[eI].emplace_back(-meshEI.first - 1, meshEJ.second, -1, add_e);
                }
                break;
            }
            case 2: {
                d_PE(mesh.vertexes[meshEI.first], mesh.vertexes[meshEJ.first], mesh.vertexes[meshEJ.second], d);
                if (d < mesh.Hhat) {
                    eePairsPerEdge[eI].emplace_back(-meshEI.first - 1, meshEJ.first, meshEJ.second, add_e);
                }
                break;
            }
            case 3: {
                d_PP(mesh.vertexes[meshEI.second], mesh.vertexes[meshEJ.first], d);
                if (d < mesh.Hhat) {
                    eePairsPerEdge[eI].emplace_back(-meshEI.second - 1, meshEJ.first, -1, add_e);
                }
                break;
            }
            case 4: {
                d_PP(mesh.vertexes[meshEI.second], mesh.vertexes[meshEJ.second], d);
                if (d < mesh.Hhat) {
                    eePairsPerEdge[eI].emplace_back(-meshEI.second - 1, meshEJ.second, -1, add_e);
                }
                break;
            }
            case 5: {
                d_PE(mesh.vertexes[meshEI.second], mesh.vertexes[meshEJ.first], mesh.vertexes[meshEJ.second], d);
                if (d < mesh.Hhat) {
                    eePairsPerEdge[eI].emplace_back(-meshEI.second - 1, meshEJ.first, meshEJ.second, add_e);
                }
                break;
            }
            case 6: {
                d_PE(mesh.vertexes[meshEJ.first], mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second], d);
                if (d < mesh.Hhat) {
                    eePairsPerEdge[eI].emplace_back(-meshEJ.first - 1, meshEI.first, meshEI.second, add_e);
                }
                break;
            }
            case 7: {
                d_PE(mesh.vertexes[meshEJ.second], mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second], d);
                if (d < mesh.Hhat) {
                    eePairsPerEdge[eI].emplace_back(-meshEJ.second - 1, meshEI.first, meshEI.second, add_e);
                }
                break;
            }
            case 8: {
                d_EE(mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second], mesh.vertexes[meshEJ.first], mesh.vertexes[meshEJ.second], d);
                if (d < mesh.Hhat) {
                    if (add_e <= -2) {
                        eePairsPerEdge[eI].emplace_back(meshEI.first, meshEI.second, meshEJ.first, -meshEJ.second - (int)mesh.surfEdges.size() - 2);
                    }
                    else {
                        eePairsPerEdge[eI].emplace_back(meshEI.first, meshEI.second, meshEJ.first, meshEJ.second);
                    }
                }
                break;
            }
            default:
                break;
            }
        }
    });

    size_t totalCCD = 0;
    for (size_t eI = 0; eI < mesh.surfEdges.size(); ++eI) {
        totalCCD += ccdPairs[eI].size();
    }

    mesh.Self_CCD_ActiveSet.reserve(mesh.Self_CCD_ActiveSet.size() + totalCCD);
    for (size_t eI = 0; eI < mesh.surfEdges.size(); ++eI) {
        mesh.Self_CCD_ActiveSet.insert(mesh.Self_CCD_ActiveSet.end(), ccdPairs[eI].begin(), ccdPairs[eI].end());
    }
}

void SpatialHash::deduplicatePPPEPairs(mesh3D& mesh,
    const std::vector<MMCVID>& ptPairs,
    const std::vector<std::vector<MMCVID>>& eePairsPerEdge)
{
    mesh.Self_EE_ActiveSet.clear();
    mesh.Self_EEeIe_ActiveSet.clear();
    std::map<MMCVID, int> constraintCounter;
    for (const auto& cI : ptPairs) {
        if (cI.data[3] < 0) {
            // PP or PE
            ++constraintCounter[cI];
        }
        else {
            mesh.Self_ActiveSet.emplace_back(cI);
        }
    }

    mesh.Self_EE_ActiveSet.clear();
    mesh.Self_EEeIe_ActiveSet.clear();
    int eI = 0;
    for (const auto& csI : eePairsPerEdge) {
        for (const auto& cI : csI) {
            if (cI.data[3] >= 0) {
                // regular EE
                mesh.Self_ActiveSet.emplace_back(cI);
            }
            else if (cI.data[3] == -1) {
                // regular PP or PE
                ++constraintCounter[cI];
            }
            else if (cI.data[3] >= -(int)mesh.surfEdges.size() - 1) {
                // nearly parallel PP or PE
                mesh.Self_EE_ActiveSet.emplace_back(cI.data[0], cI.data[1], cI.data[2], -1);
                mesh.Self_EEeIe_ActiveSet.emplace_back(eI, -cI.data[3] - 2);
            }
            else {
                // nearly parallel EE
                mesh.Self_EE_ActiveSet.emplace_back(cI.data[0], cI.data[1], cI.data[2], -cI.data[3] - (int)mesh.surfEdges.size() - 2);
                mesh.Self_EEeIe_ActiveSet.emplace_back(-1, -1);
            }
        }
        ++eI;
    }

    mesh.Self_ActiveSet.reserve(mesh.Self_ActiveSet.size() + constraintCounter.size());
    for (const auto& ccI : constraintCounter) {
        mesh.Self_ActiveSet.emplace_back(ccI.first.data[0], ccI.first.data[1], ccI.first.data[2], -ccI.second);
    }
}

void SpatialHash::calculateActivateSet(mesh3D& mesh) {
    mesh.Self_ActiveSet.clear();
    mesh.Self_EE_ActiveSet.clear();
    mesh.Self_EEeIe_ActiveSet.clear();
    mesh.Self_CCD_ActiveSet.clear();

    std::vector<MMCVID> ptPairs;
    std::vector<std::vector<MMCVID>> eePairsPerEdge;

    collectPointTrianglePairs(mesh, ptPairs);
    collectEdgeEdgePairs(mesh, eePairsPerEdge);

    deduplicatePPPEPairs(mesh, ptPairs, eePairsPerEdge);

}
