#include "collision/BVHBroadPhase.h"
#include "core/IPCDistance.h"
#include "core/SimulationParameters.h"
#include "math/Constants.h"
#include <algorithm>
#include <iostream>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

using namespace std;
using namespace Eigen;

namespace {

inline cipc::AABB bboxFromPoints(const std::vector<Eigen::Vector3d>& points)
{
    cipc::AABB box;
    for (const auto& p : points) {
        box.extend(p);
    }
    return box;
}

} // namespace

cipc::AABB BVHBroadPhase::triangleAABB(const mesh3D& mesh, int sfI)
{
    const Eigen::Vector4i& sf = mesh.surface[sfI];
    cipc::AABB box;
    box.extend(mesh.vertexes[sf[0]]);
    box.extend(mesh.vertexes[sf[1]]);
    box.extend(mesh.vertexes[sf[2]]);
    return box;
}

cipc::AABB BVHBroadPhase::triangleSweptAABB(const mesh3D& mesh, int sfI,
    const std::vector<Eigen::Vector3d>& V,
    const std::vector<Eigen::Vector3d>& searchDir, double stepSize)
{
    const Eigen::Vector4i& sf = mesh.surface[sfI];
    cipc::AABB box;
    for (int j = 0; j < 3; ++j) {
        const int vI = sf[j];
        box.extend(V[vI]);
        box.extend(V[vI] - stepSize * searchDir[vI]);
    }
    return box;
}

cipc::AABB BVHBroadPhase::edgeAABB(const mesh3D& mesh, int seI)
{
    const auto& se = mesh.surfEdges[seI];
    cipc::AABB box;
    box.extend(mesh.vertexes[se.first]);
    box.extend(mesh.vertexes[se.second]);
    return box;
}

cipc::AABB BVHBroadPhase::edgeSweptAABB(const mesh3D& mesh, int seI,
    const std::vector<Eigen::Vector3d>& V,
    const std::vector<Eigen::Vector3d>& searchDir, double stepSize)
{
    const auto& se = mesh.surfEdges[seI];
    cipc::AABB box;
    box.extend(V[se.first]);
    box.extend(V[se.first] - stepSize * searchDir[se.first]);
    box.extend(V[se.second]);
    box.extend(V[se.second] - stepSize * searchDir[se.second]);
    return box;
}

BVHBroadPhase::BVHBroadPhase(const mesh3D& mesh, double voxelSize, bool use_V_prev)
{
    build(mesh, voxelSize);
}

BVHBroadPhase::BVHBroadPhase(const mesh3D& mesh, const Eigen::VectorXd& searchDir,
    double curMaxStepSize, double voxelSize, bool use_V_prev)
{
    std::vector<Eigen::Vector3d> dir(mesh.vertexNum);
    for (int i = 0; i < mesh.vertexNum; ++i) {
        dir[i] = searchDir.segment<3>(i * 3);
    }
    build(mesh, dir, curMaxStepSize, voxelSize, use_V_prev);
}

void BVHBroadPhase::build(const mesh3D& mesh, double voxelSize)
{
    (void)voxelSize;
    isMotionBuild_ = false;
    motionStepSize_ = 0.0;
    motionSearchDir_.clear();

    triAABBs_.resize(mesh.surface.size());
    edgeAABBs_.resize(mesh.surfEdges.size());
    pointAABBs_.resize(mesh.surfVerts.size());

    tbb::parallel_for(0, (int)mesh.surface.size(), 1, [&](int sfI) {
        triAABBs_[sfI] = triangleAABB(mesh, sfI);
    });
    tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int seI) {
        edgeAABBs_[seI] = edgeAABB(mesh, seI);
    });
    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI) {
        const int vI = mesh.surfVerts[svI];
        pointAABBs_[svI] = cipc::AABB(mesh.vertexes[vI].array() - 1e-9, mesh.vertexes[vI].array() + 1e-9);
    });

    triBVH_.build(triAABBs_);
    edgeBVH_.build(edgeAABBs_);
}

void BVHBroadPhase::build(const mesh3D& mesh, const std::vector<Eigen::Vector3d>& searchDir,
    double& curMaxStepSize, double voxelSize, bool use_V_prev)
{
    (void)voxelSize;
    isMotionBuild_ = true;
    motionStepSize_ = curMaxStepSize;
    motionSearchDir_ = searchDir;

    const std::vector<Eigen::Vector3d>& V = use_V_prev ? mesh.V_prev : mesh.vertexes;

    triAABBs_.resize(mesh.surface.size());
    edgeAABBs_.resize(mesh.surfEdges.size());
    pointAABBs_.resize(mesh.surfVerts.size());

    tbb::parallel_for(0, (int)mesh.surface.size(), 1, [&](int sfI) {
        triAABBs_[sfI] = triangleSweptAABB(mesh, sfI, V, searchDir, curMaxStepSize);
    });
    tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int seI) {
        edgeAABBs_[seI] = edgeSweptAABB(mesh, seI, V, searchDir, curMaxStepSize);
    });
    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI) {
        const int vI = mesh.surfVerts[svI];
        cipc::AABB box;
        box.extend(V[vI]);
        box.extend(V[vI] - curMaxStepSize * searchDir[vI]);
        pointAABBs_[svI] = box;
    });

    triBVH_.build(triAABBs_);
    edgeBVH_.build(edgeAABBs_);
}

void BVHBroadPhase::calculateActivateSet(mesh3D& mesh)
{
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


void BVHBroadPhase::collectPointTrianglePairs(mesh3D& mesh, std::vector<MMCVID>& ptPairs)
{
    const double sqrtDHat = std::sqrt(mesh.Hhat);
    std::vector<std::vector<MMCVID>> constraintSetPT(mesh.surfVerts.size());
    std::vector<std::vector<std::pair<int, int>>> ccdPairs(mesh.surfVerts.size());

    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int svI) {
        const int vI = mesh.surfVerts[svI];
        thread_local std::vector<int> triInds;
        triInds.clear();

        cipc::AABB queryBox;
        queryBox.extend(mesh.vertexes[vI]);
        queryBox.expand(sqrtDHat);

        triBVH_.queryAABB(queryBox, [&](int sfI) {
            if (queryBox.overlaps(triAABBs_[sfI])) {
                triInds.emplace_back(sfI);
            }
        });

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

void BVHBroadPhase::collectEdgeEdgePairs(mesh3D& mesh, std::vector<std::vector<MMCVID>>& eePairsPerEdge)
{
    const double sqrtDHat = std::sqrt(mesh.Hhat);
    eePairsPerEdge.assign(mesh.surfEdges.size(), std::vector<MMCVID>());
    std::vector<std::vector<std::pair<int, int>>> ccdPairs(mesh.surfEdges.size());

    tbb::parallel_for(0, (int)mesh.surfEdges.size(), 1, [&](int eI) {
        const auto& meshEI = mesh.surfEdges[eI];
        thread_local std::vector<int> edgeInds;
        edgeInds.clear();

        cipc::AABB queryBox;
        queryBox.extend(mesh.vertexes[meshEI.first]);
        queryBox.extend(mesh.vertexes[meshEI.second]);
        queryBox.expand(sqrtDHat);

        edgeBVH_.queryAABB(queryBox, [&](int seJ) {
            if (queryBox.overlaps(edgeAABBs_[seJ]) && seJ > eI) {
                edgeInds.emplace_back(seJ);
            }
        });

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

void BVHBroadPhase::deduplicatePPPEPairs(mesh3D& mesh,
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


// Voxel-specific helpers are no-ops for BVH.
void BVHBroadPhase::locateVoxelAxisIndex(const Eigen::Vector3d& pos, Eigen::Array<int, 1, 3>& voxelAxisIndex) const
{
    (void)pos;
    voxelAxisIndex.setZero();
}

int BVHBroadPhase::locateVoxelIndex(const Eigen::Vector3d& pos) const
{
    (void)pos;
    return 0;
}

int BVHBroadPhase::voxelAxisIndex2VoxelIndex(const int voxelAxisIndex[3]) const
{
    (void)voxelAxisIndex;
    return 0;
}

int BVHBroadPhase::voxelAxisIndex2VoxelIndex(int ix, int iy, int iz) const
{
    (void)ix;
    (void)iy;
    (void)iz;
    return 0;
}

void BVHBroadPhase::queryPointForTriangles(const Eigen::Vector3d& pos, double radius, std::unordered_set<int>& triInds) const
{
    triInds.clear();
    cipc::AABB queryBox(pos.array() - radius, pos.array() + radius);
    triBVH_.queryAABB(queryBox, [&](int sfI) {
        if (queryBox.overlaps(triAABBs_[sfI])) {
            triInds.insert(sfI);
        }
    });
}

void BVHBroadPhase::queryPointForTriangles(int svI, std::unordered_set<int>& sTriInds) const
{
    sTriInds.clear();
    if (svI < 0 || svI >= (int)pointAABBs_.size()) {
        return;
    }
    cipc::AABB queryBox = pointAABBs_[svI];
    // Expand to mimic SpatialHash's swept-voxel query.
    // averageEdgeLength is not stored here; the caller already builds pointAABBs_
    // without this expansion, but the CCD step size is usually small enough that
    // the exact swept AABB is a tighter and more correct query region.
    // To strictly match SpatialHash, expand by the mesh's average edge length.
    // We leave this as the tight swept AABB for now; if numerical matching is
    // required, the motion build should expand pointAABBs_ accordingly.
    triBVH_.queryAABB(queryBox, [&](int sfI) {
        if (queryBox.overlaps(triAABBs_[sfI])) {
            sTriInds.insert(sfI);
        }
    });
}

void BVHBroadPhase::queryTriangleForPoints(const Eigen::Vector3d& v0,
    const Eigen::Vector3d& v1,
    const Eigen::Vector3d& v2,
    double radius, std::unordered_set<int>& pointInds) const
{
    pointInds.clear();
    Eigen::Vector3d leftBottom = v0.array().min(v1.array()).min(v2.array());
    Eigen::Vector3d rightTop = v0.array().max(v1.array()).max(v2.array());
    cipc::AABB queryBox(leftBottom.array() - radius, rightTop.array() + radius);
    // Surface point BVH is not maintained; return empty.
    (void)queryBox;
}

void BVHBroadPhase::queryTriangleForEdges(const Eigen::Vector3d& v0,
    const Eigen::Vector3d& v1,
    const Eigen::Vector3d& v2,
    double radius, std::unordered_set<int>& edgeInds) const
{
    edgeInds.clear();
    Eigen::Vector3d leftBottom = v0.array().min(v1.array()).min(v2.array());
    Eigen::Vector3d rightTop = v0.array().max(v1.array()).max(v2.array());
    cipc::AABB queryBox(leftBottom.array() - radius, rightTop.array() + radius);
    edgeBVH_.queryAABB(queryBox, [&](int seI) {
        if (queryBox.overlaps(edgeAABBs_[seI])) {
            edgeInds.insert(seI);
        }
    });
}

void BVHBroadPhase::queryPointForPrimitives(const Eigen::Vector3d& pos,
    const Eigen::Vector3d& dir, std::unordered_set<int>& sVInds,
    std::unordered_set<int>& sEdgeInds, std::unordered_set<int>& sTriInds) const
{
    sVInds.clear();
    sEdgeInds.clear();
    sTriInds.clear();
    (void)pos;
    (void)dir;
}

void BVHBroadPhase::queryEdgeForPE(const Eigen::Vector3d& vBegin, const Eigen::Vector3d& vEnd,
    std::vector<int>& svInds, std::vector<int>& edgeInds) const
{
    svInds.resize(0);
    edgeInds.resize(0);
    (void)vBegin;
    (void)vEnd;
}

void BVHBroadPhase::queryEdgeForEdgesWithBBoxCheck(
    const mesh3D& mesh,
    const Eigen::Matrix<double, 1, 3>& vBegin,
    const Eigen::Matrix<double, 1, 3>& vEnd,
    double radius, std::vector<int>& edgeInds,
    int eIq) const
{
    Eigen::Matrix<double, 1, 3> leftBottom = vBegin.array().min(vEnd.array()) - radius;
    Eigen::Matrix<double, 1, 3> rightTop = vBegin.array().max(vEnd.array()) + radius;
    cipc::AABB queryBox(leftBottom, rightTop);

    edgeInds.resize(0);
    edgeBVH_.queryAABB(queryBox, [&](int seJ) {
        if (seJ > eIq && queryBox.overlaps(edgeAABBs_[seJ])) {
            const auto& meshEJ = mesh.surfEdges[seJ];
            Eigen::Array<double, 1, 3> bboxEJTopRight = mesh.vertexes[meshEJ.first].array().max(mesh.vertexes[meshEJ.second].array());
            Eigen::Array<double, 1, 3> bboxEJBottomLeft = mesh.vertexes[meshEJ.first].array().min(mesh.vertexes[meshEJ.second].array());
            if (!((bboxEJBottomLeft - rightTop.array() > 0.0).any() || (leftBottom.array() - bboxEJTopRight > 0.0).any())) {
                edgeInds.emplace_back(seJ);
            }
        }
    });
    std::sort(edgeInds.begin(), edgeInds.end());
    edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
}

void BVHBroadPhase::queryEdgeForEdgesWithBBoxCheck(const mesh3D& mesh,
    const std::vector<Eigen::Vector3d>& searchDir, double curMaxStepSize,
    int seI, std::unordered_set<int>& sEdgeInds) const
{
    sEdgeInds.clear();
    if (seI < 0 || seI >= (int)mesh.surfEdges.size()) {
        return;
    }

    const Eigen::Matrix<double, 1, 3>& eI_v0 = mesh.vertexes[mesh.surfEdges[seI].first];
    const Eigen::Matrix<double, 1, 3>& eI_v1 = mesh.vertexes[mesh.surfEdges[seI].second];
    Eigen::Matrix<double, 1, 3> eI_v0t = eI_v0 - curMaxStepSize * searchDir[mesh.surfEdges[seI].first].transpose();
    Eigen::Matrix<double, 1, 3> eI_v1t = eI_v1 - curMaxStepSize * searchDir[mesh.surfEdges[seI].second].transpose();
    Eigen::Array<double, 1, 3> bboxEITopRight = eI_v0.array().max(eI_v0t.array()).max(eI_v1.array()).max(eI_v1t.array());
    Eigen::Array<double, 1, 3> bboxEIBottomLeft = eI_v0.array().min(eI_v0t.array()).min(eI_v1.array()).min(eI_v1t.array());

    cipc::AABB queryBox{ Eigen::Vector3d(bboxEIBottomLeft), Eigen::Vector3d(bboxEITopRight) };
    edgeBVH_.queryAABB(queryBox, [&](int seJ) {
        if (seJ <= seI || !queryBox.overlaps(edgeAABBs_[seJ])) {
            return;
        }
        const auto& meshEJ = mesh.surfEdges[seJ];
        Eigen::Matrix<double, 1, 3> eJ_v0 = mesh.vertexes[meshEJ.first];
        Eigen::Matrix<double, 1, 3> eJ_v1 = mesh.vertexes[meshEJ.second];
        Eigen::Matrix<double, 1, 3> eJ_v0t = eJ_v0 - curMaxStepSize * searchDir[meshEJ.first].transpose();
        Eigen::Matrix<double, 1, 3> eJ_v1t = eJ_v1 - curMaxStepSize * searchDir[meshEJ.second].transpose();
        Eigen::Array<double, 1, 3> bboxEJTopRight = eJ_v0.array().max(eJ_v0t.array()).max(eJ_v1.array()).max(eJ_v1t.array());
        Eigen::Array<double, 1, 3> bboxEJBottomLeft = eJ_v0.array().min(eJ_v0t.array()).min(eJ_v1.array()).min(eJ_v1t.array());
        if (!((bboxEJBottomLeft - bboxEITopRight > 0.0).any() || (bboxEIBottomLeft - bboxEJTopRight > 0.0).any())) {
            sEdgeInds.insert(seJ);
        }
    });
}
