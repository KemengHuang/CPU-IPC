#include "IPCSolver.h"
#include "ContactMechanics.h"
#include "FeasibleStep.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>
#include "StageTimer.h"
#include <chrono>
#include "Friction.h"
#include "NewtonLinearSystem.h"
#include "RuntimePaths.h"
#include <limits>
#include <stdexcept>


namespace {

struct NewtonWorkspace {
    NewtonWorkspace(int vertexCount, NewtonLinearSystem& sharedLinearSystem)
        : gradient(static_cast<size_t>(vertexCount), Vector3d::Zero()),
          vertexMutex(static_cast<size_t>(vertexCount)),
          linearSystem(sharedLinearSystem)
    {
    }

    vector<Vector3d> gradient;
    BHessian hessian;
    vector<tbb::spin_mutex> vertexMutex;
    NewtonLinearSystem& linearSystem;
};

struct EnergyWorkspace {
    Eigen::VectorXd constraintValues;
    Eigen::VectorXd barrierValues;
    Friction::EnergyWorkspace friction;
};

void buildConstraintStartIndsWithMM(const vector<int>& activeSet,
    const std::vector<EncodedContact>& MMActiveSet,
    std::vector<int>& constraintStartInds) {
    constraintStartInds.resize(1);
    constraintStartInds[0] = 0;
    constraintStartInds.emplace_back(
        constraintStartInds.back() + static_cast<int>(activeSet.size()));
    constraintStartInds.emplace_back(
        constraintStartInds.back() + static_cast<int>(MMActiveSet.size()));
}

} // namespace

void updateInertialTarget(mesh3D& mesh) {
    mesh.inertialTarget.resize(mesh.vertexes.size());
	double gravity = -9.8;
    if(!mesh.apply_gravity) {
        gravity = 0;
	}
    Vector3d gravityDtSq = Vector3d(0, gravity, 0) * mesh.IPC_dt * mesh.IPC_dt;

    tbb::parallel_for(0, (int)mesh.vertexes.size(), 1, [&](int vI) {
        mesh.inertialTarget[vI] =
            mesh.V_prev[vI] + mesh.velocities[vI] * mesh.IPC_dt + gravityDtSq;
        }

    );

}

namespace {

void computeEGradient(const mesh3D& mesh, vector<Vector3d>& gradient) {
    //Mass part
    Matrix3d massM;
    massM.setIdentity();

    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int vI) {
        gradient[vI] += mesh.masses[vI] * massM
            * (mesh.vertexes[vI] - mesh.inertialTarget[vI]);
        }

    );

    // Internal elastic contribution.
    vector<tbb::spin_mutex> countMutex(mesh.vertexNum);
    tbb::parallel_for(0, mesh.tetrahedraNum, 1, [&](int ii) {
        const TetPFPX& PFPX = mesh.tetrahedraPFPX[ii];
        Matrix3d F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_tetrahedra_inverse[ii];

        Matrix3d PEPF = computePEPF_StableNHK3D_2_double(F, mesh.lengthRate, mesh.volumeRate);
        Eigen::Map<const Eigen::Matrix<double, 9, 1>> pepf(PEPF.data());
        Eigen::Matrix<double, 12, 1> f;
        f.noalias() = mesh.volum[ii] * PFPX.transpose() * pepf;

        for (int i = 0; i < 4; i++) {
            countMutex[mesh.tetrahedras[ii][i]].lock();
            gradient[mesh.tetrahedras[ii][i]] +=
                mesh.IPC_dt * mesh.IPC_dt * f.template segment<3>(i * 3);
            countMutex[mesh.tetrahedras[ii][i]].unlock();
        }
        }

    );
}

int computeGradientAndHessian(
    mesh3D& mesh,
    vector<Vector3d>& gradient,
    BHessian& BH,
    const Ground& grd,
    vector<tbb::spin_mutex>& countMutex) {
    //calculate inertial gradient
    Matrix3d massM;
    massM.setIdentity();
    int collisionNum = 0;
    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int vI) {
        gradient[vI] = mesh.masses[vI] * massM
            * (mesh.vertexes[vI] - mesh.inertialTarget[vI]);
        }

    );
    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int vI) {
        gradient[vI] += (mesh.drag_coeff * mesh.masses[vI] * massM * (mesh.vertexes[vI] - mesh.V_prev[vI]));
        }

    );

#if defined USE_QUADRATIC_BENDING
    BH.H12x12.resize(mesh.tetrahedraNum + mesh.quadBendingInfo.size());
#else
    BH.H12x12.resize(mesh.tetrahedraNum + mesh.hingeBendingInfo.size());
#endif

    //calculate FEM gradient and hessian for tetrahedra mesh
    tbb::parallel_for(0, mesh.tetrahedraNum, 1, [&](int ii) {
        const TetPFPX& PFPX = mesh.tetrahedraPFPX[ii];
        Matrix3d F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_tetrahedra_inverse[ii];
        Matrix3d PEPF = computePEPF_StableNHK3D_2_double(F, mesh.lengthRate, mesh.volumeRate);

        Eigen::Map<const Eigen::Matrix<double, 9, 1>> pepf(PEPF.data());
        Eigen::Matrix<double, 12, 1> f;
        f.noalias() = mesh.volum[ii] * PFPX.transpose() * pepf;

        for (int i = 0; i < 4; i++) {
            countMutex[mesh.tetrahedras[ii][i]].lock();
            gradient[mesh.tetrahedras[ii][i]] +=
                mesh.IPC_dt * mesh.IPC_dt * f.template segment<3>(i * 3);
            countMutex[mesh.tetrahedras[ii][i]].unlock();

        }

        const Eigen::Matrix<double, 9, 9> Hq =
            project_StabbleNHK_2_H_3D(F, mesh.lengthRate, mesh.volumeRate);
        BH.H12x12[ii].noalias() = mesh.volum[ii] * mesh.IPC_dt * mesh.IPC_dt
            * PFPX.transpose() * Hq * PFPX;

        }

    );
    BH.D4Index = mesh.tetrahedras;
    //calculate FEM gradient and hessian for triangle mesh
    BH.H9x9.resize(mesh.triangleNum);
    tbb::parallel_for(0, mesh.triangleNum, 1, [&](int ii) {
        const TrianglePFPX& PFPX = mesh.trianglePFPX[ii];
        Matrix<double, 3, 2> F =
            calculateDs32D_double(mesh.vertexes, mesh.triangles[ii]) * mesh.DM_triangle_inverse[ii];
        Vector2d anisotropic_a = Vector2d(1, 0), anisotropic_b = Vector2d(0, 1);
        Matrix<double, 3, 2> PEPF = computePEPF_baraffwitkin_double(F, anisotropic_a, anisotropic_b, mesh.stretchStiffness,
            mesh.shearStiffness, mesh.strainRate);

        Eigen::Map<const Eigen::Matrix<double, 6, 1>> pepf(PEPF.data());
        Eigen::Matrix<double, 9, 1> f;
        f.noalias() = mesh.IPC_dt * mesh.IPC_dt * mesh.areas[ii]
            * PFPX.transpose() * pepf;

        for (int i = 0; i < 3; i++) {
            countMutex[mesh.triangles[ii][i]].lock();
            gradient[mesh.triangles[ii][i]] += f.template segment<3>(i * 3);
            countMutex[mesh.triangles[ii][i]].unlock();
        }

        const Eigen::Matrix<double, 6, 6> Hq = project_baraffwitkint_H_3D(F, anisotropic_a, anisotropic_b, mesh.stretchStiffness,
            mesh.shearStiffness, mesh.strainRate);
        BH.H9x9[ii].noalias() = mesh.areas[ii] * mesh.IPC_dt * mesh.IPC_dt
            * PFPX.transpose() * Hq * PFPX;

        }

    );
    BH.D3Index = mesh.triangles;
    int offset = mesh.tetrahedraNum;

#if defined USE_QUADRATIC_BENDING
    size_t bendingInfoSize = mesh.quadBendingInfo.size();
    std::vector<Vector3d>& verts = mesh.vertexes;

    vector<Vector4i> bendIndexes(bendingInfoSize);
    const double bendingScale = mesh.plateRigidity * mesh.IPC_dt * mesh.IPC_dt;

    // gradient
    tbb::parallel_for(
        0, static_cast<int>(bendingInfoSize), 1, [&](int i)
        {
            const QuadBendingInfo& info = mesh.quadBendingInfo[i];
            Eigen::Matrix<double, 12, 1> localPositions;
            for (int vertex = 0; vertex < 4; ++vertex) {
                localPositions.template segment<3>(3 * vertex) = verts[info.verts[vertex]];
            }
            const Eigen::Matrix<double, 12, 1> localGradient =
                bendingScale * info.hessianBase * localPositions;
            for (int vertex = 0; vertex < 4; ++vertex) {
                countMutex[info.verts[vertex]].lock();
                gradient[info.verts[vertex]] += localGradient.template segment<3>(3 * vertex);
                countMutex[info.verts[vertex]].unlock();
            }
        }
    );

    // hessian
    size_t H12Offset = mesh.tetrahedraNum;

    tbb::parallel_for(
        0, static_cast<int>(bendingInfoSize), 1, [&](int i)
        {
            const QuadBendingInfo& info = mesh.quadBendingInfo[i];
            Vector4i localIndex = Vector4i(info.verts[0], info.verts[1], info.verts[2], info.verts[3]);

            BH.H12x12[H12Offset + i] = bendingScale * info.hessianBase;
            bendIndexes[i] = localIndex;
        }
    );

    BH.D4Index.insert(BH.D4Index.end(), bendIndexes.begin(), bendIndexes.end());
#else
    vector<Vector4i> bendIndexes(mesh.hingeBendingInfo.size());
    const int H12Offset = mesh.tetrahedraNum;
    tbb::parallel_for(
        0, static_cast<int>(mesh.hingeBendingInfo.size()), 1, [&](int hingeIndex) {
            const HingeBendingInfo& hinge = mesh.hingeBendingInfo[hingeIndex];
            Eigen::Matrix<double, 12, 1> localGradient;
            Eigen::Matrix<double, 12, 12> localHessian;
            hingeBendingGradientAndHessian(
                hinge,
                mesh.vertexes,
                mesh.plateRigidity,
                localGradient,
                localHessian);
            IglUtils::makePD<double, 12>(localHessian);

            localGradient *= mesh.IPC_dt * mesh.IPC_dt;
            localHessian *= mesh.IPC_dt * mesh.IPC_dt;
            for (int corner = 0; corner < 4; ++corner) {
                const int vertex = hinge.vertices[corner];
                countMutex[vertex].lock();
                gradient[vertex] += localGradient.template segment<3>(3 * corner);
                countMutex[vertex].unlock();
            }
            BH.H12x12[H12Offset + hingeIndex] = localHessian;
            bendIndexes[hingeIndex] = hinge.vertices;
        });

    BH.D4Index.insert(BH.D4Index.end(), bendIndexes.begin(), bendIndexes.end());
#endif


    //calculate barrier gradient and hessian for the ground plane collision
    VectorXd constraintVals;
    offset = 0;
    Evaluate_GroundConstraintVals(grd, mesh, constraintVals, offset);
    Matrix3d nnT = grd.normal * grd.normal.transpose();

    tbb::spin_mutex groundMutex;//, countMutex3;
    tbb::parallel_for(0, (int)constraintVals.size(), 1, [&](int cI) {
        double g_b, H_b;
        compute_g_b(constraintVals[cI], mesh.Hhat, g_b);
        compute_H_b(constraintVals[cI], mesh.Hhat, H_b);
        int vI = mesh.Environment_ActiveSet[cI];
        double dist = grd.normal.dot(mesh.vertexes[vI]) - grd.D;
        gradient[vI] += mesh.Kappa * g_b * 2.0 * dist * grd.normal;
        double param = 4.0 * H_b * constraintVals[cI] + 2.0 * g_b;
        if (param > 0) {
            Matrix3d Hpg = mesh.Kappa * param * nnT;
            groundMutex.lock();
            BH.H3x3.push_back(Hpg);
            BH.D1Index.push_back(vI);
            groundMutex.unlock();
        }
        }
    );

    offset = constraintVals.size();
    Evaluate_SelfPTConstraintVals(mesh, constraintVals, offset);

    tbb::parallel_for(offset, (int)constraintVals.size(), 1, [&](int cI) {
            compute_g_b(constraintVals[cI], mesh.Hhat, constraintVals[cI]);
        }
    );
    collisionNum += compute_g_dpt(mesh, mesh.Self_ActiveSet, constraintVals, gradient, offset, mesh.Kappa);

    compute_g_dee(mesh, gradient, mesh.Hhat, mesh.Kappa);

    compute_H_dpt(mesh, BH, mesh.Hhat, mesh.Kappa);
    compute_H_dee(mesh, BH, mesh.Hhat, mesh.Kappa);

#ifdef USE_FRICTION
    Friction::addGradient(
        mesh, grd, gradient,
        mesh.Fhat * mesh.IPC_dt * mesh.IPC_dt,
        mesh.friction);
    Friction::addHessian(
        mesh, grd, BH, mesh.Fhat * mesh.IPC_dt * mesh.IPC_dt,
        mesh.friction);
#endif
    return collisionNum;
}


void computeEnergy(
    const mesh3D& mesh,
    double& energyVal,
    const Ground& gd,
    double Kappa,
    EnergyWorkspace& workspace) {
    energyVal = 0;
    energyVal += getObjEnergy_StableNHK2_3D(mesh.vertexes, mesh, mesh.lengthRate, mesh.volumeRate);
    Vector2d anisotropic_a = Vector2d(1, 0), anisotropic_b = Vector2d(0, 1);
    energyVal += getObjEnergy_baraffwitkin_3D(mesh, anisotropic_a, anisotropic_b, mesh.stretchStiffness, mesh.shearStiffness, mesh.strainRate);

#if defined USE_QUADRATIC_BENDING
    const std::vector<Vector3d>& verts = mesh.vertexes;
    double bendingEnergyVal = parallel_reduce(
        tbb::blocked_range<int>(0, mesh.quadBendingInfo.size()), 0.0,
        [&](const tbb::blocked_range<int>& rg, double temp_bendingEnergyVal) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                const QuadBendingInfo& info = mesh.quadBendingInfo[i];
                Eigen::Matrix<double, 12, 1> localPositions;
                for (int vertex = 0; vertex < 4; ++vertex) {
                    localPositions.template segment<3>(3 * vertex) = verts[info.verts[vertex]];
                }
                temp_bendingEnergyVal += localPositions.dot(info.hessianBase * localPositions);
            }
            return temp_bendingEnergyVal;
        },
        [&](double left, double right) {
            return left + right;
        });

    bendingEnergyVal *= (0.5 * mesh.plateRigidity);
    energyVal += bendingEnergyVal;
#else
    energyVal += tbb::parallel_reduce(
        tbb::blocked_range<int>(0, static_cast<int>(mesh.hingeBendingInfo.size())),
        0.0,
        [&](const tbb::blocked_range<int>& range, double localEnergy) {
            for (int hinge = range.begin(); hinge != range.end(); ++hinge) {
                localEnergy += hingeBendingEnergy(
                    mesh.hingeBendingInfo[hinge],
                    mesh.vertexes,
                    mesh.plateRigidity);
            }
            return localEnergy;
        },
        [](double left, double right) { return left + right; });
#endif

    //energyVal += getObjEnergy_AniostroI5_3D(mesh.vertexes, mesh, lengthRate, contract_ratio);
    energyVal *= mesh.IPC_dt * mesh.IPC_dt;
    double deltaE = 0;

    deltaE = parallel_reduce(
        tbb::blocked_range<int>(0, mesh.vertexNum), 0.0,
        [&](const tbb::blocked_range<int>& rg, double temp_deltaE) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                temp_deltaE += (mesh.vertexes[i] - mesh.inertialTarget[i]).squaredNorm()
                    * mesh.masses[i] / 2.0;
                temp_deltaE += ((mesh.vertexes[i] - mesh.V_prev[i]).squaredNorm() * (mesh.drag_coeff) * mesh.masses[i] / 2.0);
            }
            return temp_deltaE;
        },
        [&](double left, double right) {
            return left + right;
        }
    );

    energyVal += deltaE;
    Eigen::VectorXd& constraintVals = workspace.constraintValues;
    Eigen::VectorXd& bVals = workspace.barrierValues;
    constraintVals.resize(0);
    bVals.resize(0);
    int startCI = constraintVals.size();
    Evaluate_GroundConstraintVals(gd, mesh, constraintVals, startCI);
    bVals.conservativeResize(constraintVals.size());


    tbb::parallel_for(startCI, (int)constraintVals.size(), 1, [&](int cI) {

        compute_b(constraintVals[cI], mesh.Hhat, bVals[cI]);

        }
    );


    startCI = constraintVals.size();
    Evaluate_SelfPTConstraintVals(mesh, constraintVals, startCI);
    bVals.conservativeResize(constraintVals.size());
    tbb::parallel_for(startCI, (int)constraintVals.size(), 1, [&](int cI) {

            compute_b(constraintVals[cI], mesh.Hhat, bVals[cI]);
            int duplication = mesh.Self_ActiveSet[cI - startCI][3];
            if (duplication < -1) {
                // PP or PE, handle duplication
                bVals[cI] *= -duplication;
            }
        }

    );

    startCI = constraintVals.size();
    Evaluate_SelfEEConstraintVals(mesh, constraintVals, startCI);
    //SelfCollisionHandler<dim>::evaluateConstraints(data, paraEEMMCVIDSet.back(), constraintVals);
    bVals.conservativeResize(constraintVals.size());

    tbb::parallel_for(startCI, (int)constraintVals.size(), 1, [&](int cI) {

            const EncodedContact& MMCVIDI = mesh.Self_EE_ActiveSet[cI - startCI];
            double eps_x, e;
            if (MMCVIDI[3] >= 0) {
                // EE
                compute_eps_x(mesh, MMCVIDI[0], MMCVIDI[1], MMCVIDI[2], MMCVIDI[3], eps_x);
                compute_e(mesh.vertexes[MMCVIDI[0]], mesh.vertexes[MMCVIDI[1]], mesh.vertexes[MMCVIDI[2]],
                    mesh.vertexes[MMCVIDI[3]], eps_x, e);
            }
            else {
                // PP or PE
                const std::pair<int, int>& eIeJ = mesh.Self_EEeIe_ActiveSet[cI - startCI];
                const std::pair<int, int>& eI = mesh.surfEdges[eIeJ.first];
                const std::pair<int, int>& eJ = mesh.surfEdges[eIeJ.second];
                compute_eps_x(mesh, eI.first, eI.second, eJ.first, eJ.second, eps_x);
                compute_e(mesh.vertexes[eI.first], mesh.vertexes[eI.second], mesh.vertexes[eJ.first],
                    mesh.vertexes[eJ.second], eps_x, e);
            }
            compute_b(constraintVals[cI], mesh.Hhat, bVals[cI]);
            bVals[cI] *= e;
        }

    );
    energyVal += Kappa * bVals.sum();
#ifdef USE_FRICTION
    energyVal += Friction::energy(mesh, gd, workspace.friction);
#endif
}

void stepForward(const vector<Vector3d>& dataV0,
    const vector<Vector3d>& searchDir,
    mesh3D& mesh,
    double stepSize, bool boundary_update = false) {
    //assert(dataV0.rows() == data.V.rows());
    //assert(data.V.rows() * dim == searchDir.size());
    //assert(data.V.rows() == result.V.rows());


    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int vI) {
        if (mesh.boundaryTypes[vI] == 0 || boundary_update)
            mesh.vertexes[vI] = dataV0[vI] - stepSize * searchDir[vI];
        }

    );
}


bool checkEdgeTriIntersectionIfAny(const mesh3D& mesh,
    SpatialHash& sh) {
    Eigen::ArrayXi intersected(mesh.surface.size());
    intersected.setZero();
    tbb::parallel_for(0, (int)mesh.surface.size(), 1, [&](int sfI) {
        const Eigen::RowVector4i& sfVInd = mesh.surface[sfI].transpose();
        int coDim_sfI = 3;
        std::vector<int> sEdgeInds;
        sh.queryTriangleForEdges(mesh.vertexes[sfVInd[0]], mesh.vertexes[sfVInd[1]], mesh.vertexes[sfVInd[2]], 0.0,
            sEdgeInds);
        for (const auto& eI : sEdgeInds) {
            const auto& meshEI = mesh.surfEdges[eI];
            if (meshEI.first == sfVInd[0] || meshEI.first == sfVInd[1] || meshEI.first == sfVInd[2] ||
                meshEI.second == sfVInd[0] || meshEI.second == sfVInd[1] || meshEI.second == sfVInd[2]) {
                continue;
            }

            if (mesh.boundaryTypes[meshEI.first] >= 2 && mesh.boundaryTypes[meshEI.second] >= 2 && mesh.boundaryTypes[sfVInd[0]] >= 2 && mesh.boundaryTypes[sfVInd[1]] >= 2 && mesh.boundaryTypes[sfVInd[2]] >= 2) {
                continue;
            }

            int coDim_eI = 3;//mesh.vICoDim(meshEI.first);


            if (IglUtils::segTriIntersect(mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second],
                mesh.vertexes[sfVInd[0]], mesh.vertexes[sfVInd[1]],
                mesh.vertexes[sfVInd[2]])) {
                intersected[sfI] = 1;
                break;
            }
        }
        }

    );

    if ((intersected != 0).any()) {
        return false;
    }

    return true;
}


bool isIntersected(const Ground& grd,
    SpatialHash& sh,
    const mesh3D& mesh,
    const vector<Vector3d>& V0) {
    Eigen::VectorXd constraint_vals;
    Evaluate_GroundConstraintVals(grd, mesh, constraint_vals, 0);
    for (int vI = 0; vI < constraint_vals.size(); ++vI) {
        if (constraint_vals[vI] <= 0.0) {
            return true;
        }
    }


    if (!checkEdgeTriIntersectionIfAny(mesh, sh)) {

        return true;
    }

    return false;
}

void buildCollisionSets(mesh3D& mesh,
    SpatialHash& sh,
    const Ground& gd,
    bool rehash = true) {
    if (mesh.use_barrier) {
        if (rehash) {
            sh.build(mesh, mesh.averageEdgeLenth);
        }

        sh.calculateActivateSet(mesh);
    }
    gd.calculateActivateSet(mesh);
}

double directionInfinityNorm(const vector<Vector3d>& direction);

bool lineSearch(mesh3D& mesh,
    SpatialHash& sh,
    const Ground& gd,
    const vector<Vector3d>& searchDir,
    const vector<Vector3d>& gradient,
    double& stepSize,
    double Kappa,
    IPCStepStats& stats,
    bool diagnoseLineSearch,
    bool verbose) {
    double lastEnergyVal;
    EnergyWorkspace energyWorkspace;

    computeEnergy(mesh, lastEnergyVal, gd, Kappa, energyWorkspace);

    const double gradientDotDirection = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, mesh.vertexNum),
        0.0,
        [&](const tbb::blocked_range<int>& range, double localDot) {
            for (int vertex = range.begin(); vertex != range.end(); ++vertex) {
                localDot += gradient[vertex].dot(searchDir[vertex]);
            }
            return localDot;
        },
        [](double lhs, double rhs) { return lhs + rhs; });
    if (!std::isfinite(lastEnergyVal) || !std::isfinite(gradientDotDirection)) {
        throw std::runtime_error("line search received a non-finite energy or direction");
    }
    // The actual update is x - alpha * searchDir, so dE/dalpha at zero is
    // -gradientDotDirection. A tiny negative value is treated as zero to make
    // near-converged iterations insensitive to reduction order.
    const double predictedDecrease = (std::max)(0.0, gradientDotDirection);
    constexpr double armijoCoefficient = 0.0;

    vector<Vector3d> resultV0 = mesh.vertexes;

    if (diagnoseLineSearch) {
        constexpr double finiteDifferenceStep = 1e-4;
        auto evaluateTrialEnergy = [&](double alpha) {
            stepForward(resultV0, searchDir, mesh, alpha);
            buildCollisionSets(mesh, sh, gd);
            double energy = 0.0;
            computeEnergy(mesh, energy, gd, Kappa, energyWorkspace);
            return energy;
        };

        const double positiveEnergy = evaluateTrialEnergy(finiteDifferenceStep);
        const double negativeEnergy = evaluateTrialEnergy(-finiteDifferenceStep);
        mesh.vertexes = resultV0;
        buildCollisionSets(mesh, sh, gd);

        const double finiteDifferenceGradient =
            (positiveEnergy - negativeEnergy) / (2.0 * finiteDifferenceStep);
        const double finiteDifferenceCurvature =
            (positiveEnergy - 2.0 * lastEnergyVal + negativeEnergy)
            / (finiteDifferenceStep * finiteDifferenceStep);
        const double analyticDirectionalDerivative = -gradientDotDirection;
        const double gradientRelativeError = std::abs(
            finiteDifferenceGradient - analyticDirectionalDerivative)
            / (std::max)(1e-14, std::abs(analyticDirectionalDerivative));
        const double curvatureRatio = gradientDotDirection > 1e-14
            ? finiteDifferenceCurvature / gradientDotDirection
            : 0.0;
        std::cout << std::setprecision(17)
            << "TAYLOR h=" << finiteDifferenceStep
            << " fd1=" << finiteDifferenceGradient
            << " analytic1=" << analyticDirectionalDerivative
            << " grad_relerr=" << gradientRelativeError
            << " fd2=" << finiteDifferenceCurvature
            << " model2=" << gradientDotDirection
            << " curvature_ratio=" << curvatureRatio << std::endl;
    }

    stepForward(resultV0, searchDir, mesh, stepSize);


    int numOfIntersect = 0;
    if (mesh.use_barrier) {
        sh.build(mesh, mesh.averageEdgeLenth);
        constexpr int maxIntersectionBacktracks = 64;
        while (isIntersected(gd, sh, mesh, resultV0)) {
            if (numOfIntersect >= maxIntersectionBacktracks
                || stepSize <= std::numeric_limits<double>::epsilon()) {
                throw std::runtime_error("line search could not find a non-intersecting step");
            }
            numOfIntersect++;
            if (verbose) {
                printf("intersect\n");
            }
            stepSize /= 2.0;
            stepForward(resultV0, searchDir, mesh, stepSize);
            sh.build(mesh, mesh.averageEdgeLenth);
        }
    }

    buildCollisionSets(mesh, sh, gd, false);
    double testingE;
    computeEnergy(mesh, testingE, gd, Kappa, energyWorkspace);

    int numOfLineSearch = 0;
    const double LFStepSize = stepSize;
    auto acceptableEnergy = [&](double trialEnergy, double trialStep) {
        if (!std::isfinite(trialEnergy)) {
            return false;
        }
        const double armijoLimit = lastEnergyVal
            - armijoCoefficient * trialStep * predictedDecrease;
        return trialEnergy < armijoLimit;
    };

    constexpr int maxEnergyBacktracks = 64;
    while (!acceptableEnergy(testingE, stepSize)) {
        if (numOfLineSearch >= maxEnergyBacktracks
            || stepSize <= std::numeric_limits<double>::epsilon()
                * (std::max)(1.0, LFStepSize)) {
            const double convergenceThreshold = mesh.Newton_Solver_Threshold
                * std::sqrt(mesh.bboxDiagSize2) * mesh.IPC_dt;
            if (directionInfinityNorm(searchDir) < convergenceThreshold) {
                mesh.vertexes = resultV0;
                buildCollisionSets(mesh, sh, gd);
                stats.energyBacktracks += numOfLineSearch;
                stats.intersectionBacktracks += numOfIntersect;
                stats.lineSearchBacktracks += numOfIntersect + numOfLineSearch;
                stats.maximumEnergyBacktracksPerNewton = (std::max)(
                    stats.maximumEnergyBacktracksPerNewton, numOfLineSearch);
                if (numOfLineSearch > 2) {
                    ++stats.newtonStepsWithMoreThanTwoBacktracks;
                }
                return false;
            }
            throw std::runtime_error("line search failed to find a strict energy decrease");
        }

        if (verbose) {
            std::cout << std::setprecision(17)
                << "line-search retry=" << numOfLineSearch
                << ", alpha=" << stepSize
                << ", trialE=" << testingE
                << ", E0=" << lastEnergyVal
                << ", gTp=" << gradientDotDirection << std::endl;
        }
        stepSize *= 0.5;

        ++numOfLineSearch;

        stepForward(resultV0, searchDir, mesh, stepSize);

        buildCollisionSets(mesh, sh, gd);
        computeEnergy(mesh, testingE, gd, Kappa, energyWorkspace);
    }

    if (stepSize < LFStepSize) {
        bool needRecomputeCS = false;
        if (mesh.use_barrier) {
            constexpr int maxIntersectionBacktracks = 64;
            int safeguardIterations = 0;
            while (isIntersected(gd, sh, mesh, resultV0)) {
                if (safeguardIterations++ >= maxIntersectionBacktracks
                    || stepSize <= std::numeric_limits<double>::epsilon()) {
                    throw std::runtime_error("backtracked step remains intersecting");
                }
                stepSize /= 2.0;
                numOfIntersect++;
                if (verbose) {
                    printf("intersect\n");
                }
                stepForward(resultV0, searchDir, mesh, stepSize);
                sh.build(mesh, mesh.averageEdgeLenth);
                needRecomputeCS = true;
            }
        }
        if (needRecomputeCS) {
            buildCollisionSets(mesh, sh, gd, false);
        }
    }

    stats.energyBacktracks += numOfLineSearch;
    stats.intersectionBacktracks += numOfIntersect;
    stats.lineSearchBacktracks += numOfIntersect + numOfLineSearch;
    stats.maximumEnergyBacktracksPerNewton = (std::max)(
        stats.maximumEnergyBacktracksPerNewton, numOfLineSearch);
    if (numOfLineSearch > 2) {
        ++stats.newtonStepsWithMoreThanTwoBacktracks;
    }
    stats.acceptedStepSum += stepSize;
    stats.minimumAcceptedStep = (std::min)(stats.minimumAcceptedStep, stepSize);
    stats.maximumAcceptedStep = (std::max)(stats.maximumAcceptedStep, stepSize);

    return true;
}

void suggestKappa(double& kappa, const double& Hhat, const double& bboxDiagSize2, const double& meanMass) {
    double H_b;
    compute_H_b(1.0e-16 * bboxDiagSize2, Hhat, H_b);
    kappa = 1e13 * meanMass / (4.0e-16 * bboxDiagSize2 * H_b);
}

void upperBoundKappa(double& kappa, const double& Hhat, const double& bboxDiagSize2, const double& meanMass) {
    double H_b;
    compute_H_b(1.0e-16 * bboxDiagSize2, Hhat, H_b);
    double kappaMax = 100 * 1e13 * meanMass / (4.0e-16 * bboxDiagSize2 * H_b);
    if (kappa > kappaMax) {
        kappa = kappaMax;
    }
}

void initKappa(mesh3D& mesh, const Ground& grd, double& kappa) {
    std::vector<int> constraintStartInds;
    buildConstraintStartIndsWithMM(mesh.Environment_ActiveSet, mesh.Self_ActiveSet, constraintStartInds);

    if (constraintStartInds.back()) {
        vector<Vector3d> g_E(mesh.vertexNum, Vector3d(0, 0, 0)), g_c(mesh.vertexNum, Vector3d(0, 0, 0));

        computeEGradient(mesh, g_E);
        VectorXd constraintVal;
        int startCI = constraintStartInds[0];
        Evaluate_GroundConstraintVals(grd, mesh, constraintVal, startCI);
        for (int cI = startCI; cI < constraintStartInds[1]; ++cI) {
            compute_g_b(constraintVal[cI], mesh.Hhat, constraintVal[cI]);
            int vI = mesh.Environment_ActiveSet[cI];
            double dist = grd.normal.dot(mesh.vertexes[vI]) - grd.D;
            g_c[vI] += constraintVal[cI] * 2.0 * dist * grd.normal;
        }

        startCI = constraintStartInds[1];
        Evaluate_SelfPTConstraintVals(mesh, constraintVal, startCI);

        for (int cI = startCI; cI < constraintVal.size(); ++cI) {
            compute_g_b(constraintVal[cI], mesh.Hhat, constraintVal[cI]);
        }
        compute_g_dpt(mesh, mesh.Self_ActiveSet, constraintVal, g_c, startCI, 1);
        double gsum = 0, gsnorm = 0;
        for (int i = 0; i < mesh.vertexNum; i++) {
            gsum += g_c[i].dot(g_E[i]);
            gsnorm += g_c[i].squaredNorm();
        }
        // balance current gradient at constrained DOF
        double minKappa = -gsum / gsnorm;
        if (minKappa > 0.0) {
            kappa = minKappa;
        }
        suggestKappa(minKappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
        if (kappa < minKappa) {
            kappa = minKappa;
        }
        upperBoundKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
    }

}


void postLineSearch(mesh3D& mesh, const Ground& grd, double& kappa) {
    if (kappa == 0.0) {
        initKappa(mesh, grd, kappa);
    }
    else {
        bool updateKappa = false;
        for (int i = 0; i < mesh.closeConstraintID.size(); ++i) {

            double d = grd.calculateGapFromObj(mesh, mesh.closeConstraintID[i]);
            if (d <= mesh.closeConstraintVal[i]) {
                updateKappa = true;
                break;
            }
        }
        if (!updateKappa) {
            for (int i = 0; i < mesh.closeMConstraintID.size(); ++i) {

                double d = SelfConstraintVal(mesh, mesh.closeMConstraintID[i]);

                if (d <= mesh.closeMConstraintVal[i]) {
                    updateKappa = true;
                    break;
                }
            }
        }
        if (updateKappa) {
            kappa *= 2.0;
            upperBoundKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
        }

        mesh.closeConstraintID.resize(0);
        mesh.closeMConstraintID.resize(0);
        mesh.closeConstraintVal.resize(0);
        mesh.closeMConstraintVal.resize(0);
        Eigen::VectorXd constraintVal;

        int constraintValIndStart = constraintVal.size();
        Evaluate_GroundConstraintVals(grd, mesh, constraintVal, constraintValIndStart);

        for (int i = 0; i < mesh.Environment_ActiveSet.size(); ++i) {
            if (constraintVal[constraintValIndStart + i] < mesh.dTol) {
                mesh.closeConstraintID.emplace_back(mesh.Environment_ActiveSet[i]);
                mesh.closeConstraintVal.emplace_back(constraintVal[constraintValIndStart + i]);
            }
        }

        constraintValIndStart = constraintVal.size();
        Evaluate_SelfPTConstraintVals(mesh, constraintVal, constraintValIndStart);

        for (int i = 0; i < mesh.Self_ActiveSet.size(); ++i) {
            if (constraintVal[constraintValIndStart + i] < mesh.dTol) {
                mesh.closeMConstraintID.emplace_back(mesh.Self_ActiveSet[i]);
                mesh.closeMConstraintVal.emplace_back(constraintVal[constraintValIndStart + i]);
            }
        }
    }


}

double directionInfinityNorm(const vector<Vector3d>& direction) {
    return parallel_reduce(
        tbb::blocked_range<int>(0, static_cast<int>(direction.size())), 0.0,
        [&](const tbb::blocked_range<int>& rg, double temp_max) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                for (int jj = 0; jj < 3; jj++) {
                    if (temp_max < std::abs(direction[i][jj])) {
                        temp_max = std::abs(direction[i][jj]);
                    }
                }
            }
            return temp_max;
        },
        [&](double left, double right) {
            return left > right ? left : right;
        }
    );
}

int solveBarrierSubproblem(
    mesh3D& mesh,
    SpatialHash& sh,
    Ground& gd,
    double Kappa,
    IPCStepStats& stats,
    bool verbose,
    bool diagnoseLineSearch,
    const LinearSolverOptions& linearSolverOptions,
    NewtonLinearSystem& linearSystem) {
    constexpr int iterationLimit = 10000;
    int k = 0;

    vector<Vector3d> moveDir(mesh.vertexNum, Vector3d(0, 0, 0));
    const std::size_t analysesBefore = linearSystem.symbolicAnalysisCount();
    const std::size_t factorizationsBefore = linearSystem.numericFactorizationCount();
    const double pardisoAnalysisBefore = linearSystem.pardisoAnalysisMilliseconds();
    const double pardisoFactorizationBefore = linearSystem.pardisoFactorizationMilliseconds();
    const double pardisoSolveBefore = linearSystem.pardisoSolveMilliseconds();
    NewtonWorkspace workspace(mesh.vertexNum, linearSystem);
    double beta = 1;
	int semi_implicit_kmin = 6;
    for (; k < iterationLimit; ++k) {

        StageTimer timer0, timer1, timer2, timer3, timer4;

        vector<Vector3d>& gradient = workspace.gradient;
        BHessian& BH = workspace.hessian;
        BH.clear();
        timer0.start();
        stats.collisions += computeGradientAndHessian(
            mesh, gradient, BH, gd, workspace.vertexMutex);

        timer0.stop();
        timer1.start();

        workspace.linearSystem.solve(
            mesh, BH, gradient, moveDir,
            linearSolverOptions, stats.matrixNonZeros);

        timer1.stop();

        // Test the direction solved from the current gradient/Hessian before
        // spending work on CCD or strict-energy line search.
        const double directionNorm = directionInfinityNorm(moveDir);
        const double convergenceThreshold = mesh.Newton_Solver_Threshold
            * std::sqrt(mesh.bboxDiagSize2) * mesh.IPC_dt;
        if (directionNorm < convergenceThreshold) {
            stats.assemblyMilliseconds += timer0.milliseconds();
            stats.linearSolveMilliseconds += timer1.milliseconds();
            break;
        }

        timer2.start();

        double alpha = 1.0;
        constexpr double environmentSlackness = 0.8;
        constexpr double selfCollisionSlackness = 0.8;

        limitStepByGround(
            mesh, gd, moveDir, environmentSlackness, alpha);

        if (alpha <= 0) {
            alpha = 1;
        }
        if (mesh.use_barrier) {
            double partialCCDAlpha = alpha;
            Self_largestFeasibleStepSize(
                mesh, moveDir, selfCollisionSlackness, partialCCDAlpha);

            const double maxSurfaceMove = tbb::parallel_reduce(
                tbb::blocked_range<int>(0, static_cast<int>(mesh.surfVerts.size())),
                0.0,
                [&](const tbb::blocked_range<int>& range, double localMax) {
                    for (int i = range.begin(); i != range.end(); ++i) {
                        localMax = (std::max)(
                            localMax,
                            moveDir[static_cast<size_t>(mesh.surfVerts[i])].norm());
                    }
                    return localMax;
                },
                [](double lhs, double rhs) { return (std::max)(lhs, rhs); });
            const double alpha_CFL = maxSurfaceMove > 0.0
                ? std::sqrt(mesh.Hhat) / (maxSurfaceMove * 2.0)
                : 1.0;

            double fullCCD_alpha = alpha;
            sh.build(mesh, moveDir, fullCCD_alpha, mesh.averageEdgeLenth);
            Self_largestFeasibleStepSize_CCD(
                mesh, sh, moveDir, selfCollisionSlackness, fullCCD_alpha);

            // Preserve the original IPC CFL strategy. Self_CCD_ActiveSet and
            // full-CCD candidates are now normalized by backend-independent
            // exact AABB tests, so this branch no longer depends on whether
            // SpatialHash or LBVH produced the broad candidate superset.
            alpha = (std::min)(alpha, alpha_CFL);
            if (partialCCDAlpha > 2.0 * alpha_CFL) {
                alpha = (std::min)(partialCCDAlpha, fullCCD_alpha);
                alpha = (std::max)(alpha, alpha_CFL);
            }
        }

        timer2.stop();
        timer3.start();

        const bool acceptedStep = lineSearch(
            mesh, sh, gd, moveDir, gradient, alpha, Kappa, stats,
            diagnoseLineSearch, verbose);
        timer3.stop();
        if (!acceptedStep) {
            stats.assemblyMilliseconds += timer0.milliseconds();
            stats.linearSolveMilliseconds += timer1.milliseconds();
            stats.ccdMilliseconds += timer2.milliseconds();
            stats.lineSearchMilliseconds += timer3.milliseconds();
            break;
        }
        timer4.start();
        postLineSearch(mesh, gd, Kappa);
        timer4.stop();

        float time00, time11, time22, time33, time44;
        time00 = static_cast<float>(timer0.milliseconds());
        time11 = static_cast<float>(timer1.milliseconds());
        time22 = static_cast<float>(timer2.milliseconds());
        time33 = static_cast<float>(timer3.milliseconds());
        time44 = static_cast<float>(timer4.milliseconds());

        stats.assemblyMilliseconds += time00;
        stats.linearSolveMilliseconds += time11;
        stats.ccdMilliseconds += time22;
        stats.lineSearchMilliseconds += time33;
        stats.postLineSearchMilliseconds += time44;

        if (k + 1 >= semi_implicit_kmin)
        {
            beta = (1 - alpha) * beta;
        }
        else
        {
            beta = beta;
        }

        if (beta <= convergenceThreshold)
        {
            break;
        }

    }
    if (verbose) {
        printf("newton iteration:  %d    and    Kappa:  %f\n", k, mesh.Kappa);
    }
    stats.symbolicAnalyses += workspace.linearSystem.symbolicAnalysisCount() - analysesBefore;
    stats.numericFactorizations +=
        workspace.linearSystem.numericFactorizationCount() - factorizationsBefore;
    stats.pardisoAnalysisMilliseconds +=
        workspace.linearSystem.pardisoAnalysisMilliseconds() - pardisoAnalysisBefore;
    stats.pardisoFactorizationMilliseconds +=
        workspace.linearSystem.pardisoFactorizationMilliseconds() - pardisoFactorizationBefore;
    stats.pardisoSolveMilliseconds +=
        workspace.linearSystem.pardisoSolveMilliseconds() - pardisoSolveBefore;
    stats.linearSolverThreads = (std::max)(
        stats.linearSolverThreads, workspace.linearSystem.pardisoThreadCount());
    const int pardisoFactorNonZeros = workspace.linearSystem.pardisoFactorNonZeros();
    if (pardisoFactorNonZeros > 0) {
        stats.factorNonZeros = (std::max)(
            stats.factorNonZeros,
            static_cast<std::size_t>(pardisoFactorNonZeros));
    }
    stats.newtonIterations += k;
    return k;
}

void updateVelocity(const vector<Vector3d>& currentPos, const vector<Vector3d> originalPos, vector<Vector3d>& velocity, const double& delta_t, const int& number, bool is_quasi_static) {
    tbb::parallel_for(0, number, 1, [&](int i)
        {
            if (is_quasi_static)
                velocity[i] = Vector3d(0, 0, 0);
            else
                velocity[i] = (currentPos[i] - originalPos[i]) / delta_t;
        }
    );
}

} // namespace

int solveIPCStep(
    int& stepId,
    mesh3D& mesh,
    SpatialHash& sh,
    Ground& gd,
    IPCSolverContext& context) {
    if (!context.checkpointLoadAttempted) {
        if (mesh.resumedFromCheckpoint) {
            const std::string fileVertex = RuntimePaths::tempFile("timeCost.txt");
            ifstream ifs(fileVertex);
            if (ifs) {
                ifs >> context.cumulativeStageMilliseconds[0]
                    >> context.cumulativeStageMilliseconds[1]
                    >> context.cumulativeStageMilliseconds[2]
                    >> context.cumulativeStageMilliseconds[3]
                    >> context.cumulativeStageMilliseconds[4]
                    >> context.cumulativeStepMilliseconds
                    >> context.totalNewtonIterations
                    >> context.stepIndex
                    >> mesh.Kappa;
            }
        }
        context.checkpointLoadAttempted = true;
    }

    stepId = context.stepIndex;
    IPCStepStats stats;
    stats.frame = stepId;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    upperBoundKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
    if (mesh.Kappa < 1e-16) {
        suggestKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
    }
    initKappa(mesh, gd, mesh.Kappa);
#ifdef USE_FRICTION
    Friction::initialize(mesh, gd);
#endif

    if (mesh.update_hard_constraint_functor != nullptr) {
        vector<Vector3d> moveDir(mesh.vertexNum, Vector3d(0, 0, 0));
        double new_alpha = 1;
        int boundary_vertex_num = mesh.boundary_vertexes_indices.size();
        tbb::parallel_for(0, boundary_vertex_num, 1, [&](int i)
            {
                moveDir[mesh.boundary_vertexes_indices[i]] = mesh.update_hard_constraint_functor(
                    mesh.vertexes[mesh.boundary_vertexes_indices[i]], new_alpha, mesh.IPC_dt);
            }
        );
        sh.build(mesh, moveDir, new_alpha, mesh.averageEdgeLenth);
        Self_largestFeasibleStepSize_CCD(mesh, sh, moveDir, 0.8, new_alpha);
        tbb::parallel_for(0, boundary_vertex_num, 1, [&](int i)
            {
                moveDir[mesh.boundary_vertexes_indices[i]] = mesh.update_hard_constraint_functor(
                    mesh.vertexes[mesh.boundary_vertexes_indices[i]], new_alpha, mesh.IPC_dt);
            }
        );
        vector<Vector3d> resultV0 = mesh.vertexes;
        stepForward(resultV0, moveDir, mesh, 1, true);

        sh.build(mesh, mesh.averageEdgeLenth);
        int numOfIntersect = 0;

        while (isIntersected(gd, sh, mesh, mesh.vertexes)) {
            if (numOfIntersect++ >= 64
                || new_alpha <= std::numeric_limits<double>::epsilon()) {
                throw std::runtime_error("animated boundary update remains intersecting");
            }
            new_alpha /= 2.0;
            tbb::parallel_for(0, boundary_vertex_num, 1, [&](int i)
                {
                    moveDir[mesh.boundary_vertexes_indices[i]] = mesh.update_hard_constraint_functor(
                        mesh.vertexes[mesh.boundary_vertexes_indices[i]], new_alpha, mesh.IPC_dt);
                }
            );
            stepForward(resultV0, moveDir, mesh, 1, true);
            sh.build(mesh, mesh.averageEdgeLenth);
        }
        if (context.verbose) {
            cout << "new_alpha:" << new_alpha << endl;
        }
        buildCollisionSets(mesh, sh, gd, false);
    }
    int k = 0;
    if (!context.linearSystem
        || context.linearSystem->vertexCount() != mesh.vertexNum) {
        context.linearSystem = std::make_shared<NewtonLinearSystem>(mesh.vertexNum);
    }
    constexpr int maxKappaIterations = 64;
    for (; stats.kappaIterations < maxKappaIterations; ++stats.kappaIterations) {
        mesh.closeConstraintID.resize(0);
        mesh.closeMConstraintID.resize(0);
        mesh.closeConstraintVal.resize(0);
        mesh.closeMConstraintVal.resize(0);
        k = solveBarrierSubproblem(
            mesh, sh, gd, mesh.Kappa, stats,
            context.verbose, context.diagnoseLineSearch, context.linearSolver,
            *context.linearSystem);

        VectorXd constraintVals;
        int offset = 0;
        Evaluate_GroundConstraintVals(gd, mesh, constraintVals, offset);
        offset = constraintVals.size();
        Evaluate_SelfPTConstraintVals(mesh, constraintVals, offset);

        if (constraintVals.size()) {
            const double minimumConstraint = constraintVals.minCoeff();
            const double maximumConstraint = constraintVals.maxCoeff();
            stats.minConstraintDistance2 = minimumConstraint;
            if (minimumConstraint < mesh.dTol) {
                break;
            }
            else if (maximumConstraint < mesh.Hhat) {
                break;
            }
        }
        else {
            break;
        }
    }

    if (stats.kappaIterations == maxKappaIterations) {
        throw std::runtime_error("IPC kappa loop exceeded its iteration limit");
    }
    ++stats.kappaIterations;


    updateVelocity(mesh.vertexes, mesh.V_prev, mesh.velocities, mesh.IPC_dt, mesh.vertexNum, mesh.is_quasi_static);


    mesh.V_prev = mesh.vertexes;
    updateInertialTarget(mesh);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    stats.stepMilliseconds = std::chrono::duration<double, std::milli>(end - begin).count();
    stats.kappa = mesh.Kappa;
    stats.groundContacts = mesh.Environment_ActiveSet.size();
    stats.selfContacts = mesh.Self_ActiveSet.size();
    stats.mollifiedContacts = mesh.Self_EE_ActiveSet.size();

    context.cumulativeStageMilliseconds[0] += stats.assemblyMilliseconds;
    context.cumulativeStageMilliseconds[1] += stats.linearSolveMilliseconds;
    context.cumulativeStageMilliseconds[2] += stats.ccdMilliseconds;
    context.cumulativeStageMilliseconds[3] += stats.lineSearchMilliseconds;
    context.cumulativeStageMilliseconds[4] += stats.postLineSearchMilliseconds;
    context.cumulativeStepMilliseconds += stats.stepMilliseconds;
    context.totalNewtonIterations += stats.newtonIterations;
    context.totalCollisions += stats.collisions;
    context.lastStep = stats;

    const double averageCollision = context.totalNewtonIterations > 0
        ? context.totalCollisions / context.totalNewtonIterations
        : 0.0;

    if (context.verbose) {
        std::cout << "finished step " << stats.frame
            << " in " << stats.stepMilliseconds << " ms"
            << ", Newton=" << stats.newtonIterations
            << ", kappa loops=" << stats.kappaIterations
            << ", nnz=" << stats.matrixNonZeros
            << ", analyze=" << stats.symbolicAnalyses << std::endl;
        std::cout << "line-search energy retries=" << stats.energyBacktracks
            << " (mean/Newton="
            << (stats.newtonIterations > 0
                ? static_cast<double>(stats.energyBacktracks) / stats.newtonIterations
                : 0.0)
            << ", max/Newton=" << stats.maximumEnergyBacktracksPerNewton
            << ", Newton>2=" << stats.newtonStepsWithMoreThanTwoBacktracks
            << "), intersection retries=" << stats.intersectionBacktracks << std::endl;
        std::cout << "totalCollision=" << context.totalCollisions
            << ", averageCollision=" << averageCollision
            << ", total=" << context.cumulativeStepMilliseconds << " ms" << std::endl;
        if (stats.linearSolverThreads > 0) {
            std::cout << "PARDISO phase ms: analyze=" << stats.pardisoAnalysisMilliseconds
                << ", factorize=" << stats.pardisoFactorizationMilliseconds
                << ", solve=" << stats.pardisoSolveMilliseconds
                << ", threads=" << stats.linearSolverThreads
                << ", factor nnz=" << stats.factorNonZeros << std::endl;
        }
    }

    ++context.stepIndex;

    if (context.writeRuntimeFiles) {
        RuntimePaths::initialize();
        ofstream outTime(RuntimePaths::outputFile("timeCost.txt"));
        outTime << "time0: " << context.cumulativeStageMilliseconds[0] / 1000.0 << endl;
        outTime << "time1: " << context.cumulativeStageMilliseconds[1] / 1000.0 << endl;
        outTime << "time2: " << context.cumulativeStageMilliseconds[2] / 1000.0 << endl;
        outTime << "time3: " << context.cumulativeStageMilliseconds[3] / 1000.0 << endl;
        outTime << "time4: " << context.cumulativeStageMilliseconds[4] / 1000.0 << endl;
        outTime << "timeAll: " << context.cumulativeStepMilliseconds / 1000.0 << endl;
        outTime << "total iter: " << context.totalNewtonIterations << endl;
        outTime << "frames: " << context.stepIndex << endl;
        outTime << "totalCollisionNum: " << context.totalCollisions << endl;
        outTime << "averageCollision: " << averageCollision << endl;

        bool metricsHasContent = false;
        if (mesh.resumedFromCheckpoint) {
            ifstream existingMetrics(RuntimePaths::outputFile("metrics.csv"));
            metricsHasContent = existingMetrics && existingMetrics.peek() != ifstream::traits_type::eof();
        }
        std::ios_base::openmode metricsMode = std::ios::out | std::ios::app;
        if (!context.metricsInitialized && !mesh.resumedFromCheckpoint) {
            metricsMode = std::ios::out | std::ios::trunc;
        }
        ofstream metrics(RuntimePaths::outputFile("metrics.csv"), metricsMode);
        if (!context.metricsInitialized && !metricsHasContent) {
            metrics << "frame,step_ms,assembly_ms,linear_ms,ccd_ms,line_search_ms,post_ms,"
                "newton,kappa_loops,backtracks,energy_backtracks,intersection_backtracks,"
                "max_energy_backtracks_per_newton,newton_steps_over_two_backtracks,"
                "mean_alpha,min_alpha,max_alpha,collisions,kappa,min_distance2,"
                "ground_contacts,self_contacts,mollified_contacts,matrix_nnz,"
                "symbolic_analyses,numeric_factorizations,pardiso_analysis_ms,"
                "pardiso_factorization_ms,pardiso_solve_ms,linear_solver_threads,"
                "factor_nnz\n";
        }
        metrics << std::setprecision(17)
            << stats.frame << ',' << stats.stepMilliseconds << ','
            << stats.assemblyMilliseconds << ',' << stats.linearSolveMilliseconds << ','
            << stats.ccdMilliseconds << ',' << stats.lineSearchMilliseconds << ','
            << stats.postLineSearchMilliseconds << ',' << stats.newtonIterations << ','
            << stats.kappaIterations << ',' << stats.lineSearchBacktracks << ','
            << stats.energyBacktracks << ',' << stats.intersectionBacktracks << ','
            << stats.maximumEnergyBacktracksPerNewton << ','
            << stats.newtonStepsWithMoreThanTwoBacktracks << ','
            << (stats.newtonIterations > 0
                ? stats.acceptedStepSum / stats.newtonIterations
                : 0.0) << ','
            << stats.minimumAcceptedStep << ',' << stats.maximumAcceptedStep << ','
            << stats.collisions << ',' << stats.kappa << ',' << stats.minConstraintDistance2 << ','
            << stats.groundContacts << ',' << stats.selfContacts << ','
            << stats.mollifiedContacts << ',' << stats.matrixNonZeros << ','
            << stats.symbolicAnalyses << ',' << stats.numericFactorizations << ','
            << stats.pardisoAnalysisMilliseconds << ','
            << stats.pardisoFactorizationMilliseconds << ','
            << stats.pardisoSolveMilliseconds << ','
            << stats.linearSolverThreads << ',' << stats.factorNonZeros << '\n';
        context.metricsInitialized = true;

        if (context.writeCheckpoints && context.stepIndex % 10 == 0) {
            ofstream outTime2(RuntimePaths::tempFile("timeCost.txt"));
            outTime2 << context.cumulativeStageMilliseconds[0] << endl;
            outTime2 << context.cumulativeStageMilliseconds[1] << endl;
            outTime2 << context.cumulativeStageMilliseconds[2] << endl;
            outTime2 << context.cumulativeStageMilliseconds[3] << endl;
            outTime2 << context.cumulativeStageMilliseconds[4] << endl;
            outTime2 << context.cumulativeStepMilliseconds << endl;
            outTime2 << context.totalNewtonIterations << endl;
            outTime2 << context.stepIndex << endl;
            outTime2 << mesh.Kappa << endl;
            if (!mesh.output_tetTempData()) {
                std::cerr << "Warning: failed to save simulation checkpoint" << std::endl;
            }
        }
    }

    return stats.newtonIterations;
}
