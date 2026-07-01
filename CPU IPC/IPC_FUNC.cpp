#include "IPC_FUNC.h"
#include "IPCdistanceFuncs.h"
#include "fem_parameters.h"
#include "IPCtimeStepFuns.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>
#include "fem_timer.h"
#include <chrono>
#include "FrictionUtils.hpp"
#include "Solver.h"
#include "SimulationParameters.h"
#include <assert.h>


using namespace FEM;
using namespace std;
using namespace Eigen;
#define NEWB2
//index count to export
int step_index = 0;
double time_total = 0;

double ttime0 = 0;
double ttime1 = 0;
double ttime2 = 0;
double ttime3 = 0;
double ttime4 = 0;
double totalCollision = 0;
int total_iter = 0;

void buildConstraintStartIndsWithMM(const vector<int>& activeSet,
    const std::vector<MMCVID>& MMActiveSet,
    std::vector<int>& constraintStartInds) {
    constraintStartInds.resize(1);
    constraintStartInds[0] = 0;
    constraintStartInds.emplace_back(constraintStartInds.back() + activeSet.size());
    constraintStartInds.emplace_back(constraintStartInds.back() + MMActiveSet.size());
}

void buildFriction(mesh3D& mesh, const Ground& grd) {
    Eigen::VectorXd constraintVals;

    {
        int startCI = 0;
        Evaluate_SelfPTConstraintVals(mesh, constraintVals, 0);
        mesh.Self_lambda_lastH.resize(constraintVals.size() - startCI);
        //TODO: parallelize
        for (int i = 0; i < mesh.Self_lambda_lastH.size(); ++i) {
            compute_g_b(constraintVals[startCI + i], mesh.Hhat, mesh.Self_lambda_lastH[i]);
            mesh.Self_lambda_lastH[i] *= -mesh.Kappa * 2.0 * std::sqrt(constraintVals[startCI + i]);
            if (mesh.Self_ActiveSet[i][3] < -1) {
                // PP or PE duplication
                mesh.Self_lambda_lastH[i] *= -mesh.Self_ActiveSet[i][3];
            }
        }

        mesh.MMDistCoord.resize(mesh.Self_ActiveSet.size());
        mesh.MMTanBasis.resize(mesh.Self_ActiveSet.size());
        for (int cI = 0; cI < mesh.Self_ActiveSet.size(); ++cI) {
            const auto& MMCVIDI = mesh.Self_ActiveSet[cI];
            if (MMCVIDI[0] >= 0) {
                // edge-edge
                IPC::computeClosestPoint_EE(mesh.vertexes[MMCVIDI[0]], mesh.vertexes[MMCVIDI[1]],
                    mesh.vertexes[MMCVIDI[2]], mesh.vertexes[MMCVIDI[3]], mesh.MMDistCoord[cI]);
                IPC::computeTangentBasis_EE(mesh.vertexes[MMCVIDI[0]], mesh.vertexes[MMCVIDI[1]],
                    mesh.vertexes[MMCVIDI[2]], mesh.vertexes[MMCVIDI[3]], mesh.MMTanBasis[cI]);
            }
            else {
                // point-triangle and degenerate edge-edge
                assert(MMCVIDI[1] >= 0);
                if (MMCVIDI[2] < 0) {
                    // PP
                    mesh.MMDistCoord[cI].setZero(); // Store something instead of random memory
                    IPC::computeTangentBasis_PP(mesh.vertexes[-MMCVIDI[0] - 1], mesh.vertexes[MMCVIDI[1]],
                        mesh.MMTanBasis[cI]);
                }
                else if (MMCVIDI[3] < 0) {
                    // PE
                    IPC::computeClosestPoint_PE(mesh.vertexes[-MMCVIDI[0] - 1],
                        mesh.vertexes[MMCVIDI[1]], mesh.vertexes[MMCVIDI[2]],
                        mesh.MMDistCoord[cI][0]);
                    mesh.MMDistCoord[cI][1] = 0; // Store something instead of random memory
                    IPC::computeTangentBasis_PE(mesh.vertexes[-MMCVIDI[0] - 1],
                        mesh.vertexes[MMCVIDI[1]], mesh.vertexes[MMCVIDI[2]],
                        mesh.MMTanBasis[cI]);
                }
                else {
                    // PT
                    IPC::computeClosestPoint_PT(mesh.vertexes[-MMCVIDI[0] - 1],
                        mesh.vertexes[MMCVIDI[1]], mesh.vertexes[MMCVIDI[2]],
                        mesh.vertexes[MMCVIDI[3]],
                        mesh.MMDistCoord[cI]);
                    IPC::computeTangentBasis_PT(mesh.vertexes[-MMCVIDI[0] - 1],
                        mesh.vertexes[MMCVIDI[1]], mesh.vertexes[MMCVIDI[2]],
                        mesh.vertexes[MMCVIDI[3]],
                        mesh.MMTanBasis[cI]);
                }
            }
        }
        mesh.Self_activeSet_lastH = mesh.Self_ActiveSet;
    }

    // ground
    {
        int startCI = constraintVals.size();
        Evaluate_GroundConstraintVals(grd, mesh, constraintVals, startCI);
        mesh.Environment_lambda_lastH.resize(constraintVals.size() - startCI);
        for (int i = 0; i < mesh.Environment_lambda_lastH.size(); ++i) {
            compute_g_b(constraintVals[startCI + i], mesh.Hhat, mesh.Environment_lambda_lastH[i]);
            mesh.Environment_lambda_lastH[i] *= -mesh.Kappa * 2.0 * std::sqrt(constraintVals[startCI + i]);
        }

        mesh.Environment_activeSet_lastH = mesh.Environment_ActiveSet;
    }
}


void compute_fiction_gradient(mesh3D& mesh, const Ground& grd, std::vector<Eigen::Vector3d>& grad_inc, double eps2,
    double coef) {
    double eps = std::sqrt(eps2);
    //TODO: parallelize
    for (int cI = 0; cI < mesh.Self_activeSet_lastH.size(); ++cI) {
        Eigen::RowVector3d relDX3D;

        const auto& MMCVIDI = mesh.Self_activeSet_lastH[cI];
        if (MMCVIDI[0] >= 0) {
            // edge-edge
            IPC::computeRelDX_EE(mesh.vertexes[MMCVIDI[0]] - mesh.V_prev[MMCVIDI[0]],
                mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]],
                mesh.vertexes[MMCVIDI[2]] - mesh.V_prev[MMCVIDI[2]],
                mesh.vertexes[MMCVIDI[3]] - mesh.V_prev[MMCVIDI[3]],
                mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], relDX3D);

            Eigen::Vector2d relDX = (relDX3D * mesh.MMTanBasis[cI]).transpose();
            double relDXSqNorm = relDX.squaredNorm();
            if (relDXSqNorm > eps2) {
                relDX /= std::sqrt(relDXSqNorm);
            }
            else {
                double f1_div_relDXNorm;
                IPC::f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);
                relDX *= f1_div_relDXNorm;
            }

            Eigen::Matrix<double, 12, 1> TTTDX;
            IPC::liftRelDXTanToMesh_EE(relDX, mesh.MMTanBasis[cI],
                mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], TTTDX);
            TTTDX *= coef * mesh.Self_lambda_lastH[cI];
            //            for (int i = 0; i < 4; i++) {
            //                grad_inc[MMCVIDI[i]] += TTTDX.template segment<3>(3 * i);
            //            }

            grad_inc[MMCVIDI[0]] += TTTDX.template segment<3>(0);
            grad_inc[MMCVIDI[1]] += TTTDX.template segment<3>(3);
            grad_inc[MMCVIDI[2]] += TTTDX.template segment<3>(6);
            grad_inc[MMCVIDI[3]] += TTTDX.template segment<3>(9);
            //            grad_inc.template segment<3>(MMCVIDI[0] * 3) += TTTDX.template segment<3>(0);
            //            grad_inc.template segment<3>(MMCVIDI[1] * 3) += TTTDX.template segment<3>(3);
            //            grad_inc.template segment<3>(MMCVIDI[2] * 3) += TTTDX.template segment<3>(6);
            //            grad_inc.template segment<3>(MMCVIDI[3] * 3) += TTTDX.template segment<3>(9);
        }
        else {
            // point-triangle and degenerate edge-edge
            assert(MMCVIDI[1] >= 0);
            if (MMCVIDI[2] < 0) {
                // PP
                IPC::computeRelDX_PP(mesh.vertexes[-MMCVIDI[0] - 1] - mesh.V_prev[-MMCVIDI[0] - 1],
                    mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]], relDX3D);

                Eigen::Vector2d relDX = (relDX3D * mesh.MMTanBasis[cI]).transpose();
                double relDXSqNorm = relDX.squaredNorm();
                if (relDXSqNorm > eps2) {
                    relDX /= std::sqrt(relDXSqNorm);
                }
                else {
                    double f1_div_relDXNorm;
                    IPC::f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);
                    relDX *= f1_div_relDXNorm;
                }

                Eigen::Matrix<double, 6, 1> TTTDX;
                IPC::liftRelDXTanToMesh_PP(relDX, mesh.MMTanBasis[cI], TTTDX);
                TTTDX *= coef * mesh.Self_lambda_lastH[cI];

                grad_inc[-MMCVIDI[0] - 1] += TTTDX.template segment<3>(0);
                grad_inc[MMCVIDI[1]] += TTTDX.template segment<3>(3);
            }
            else if (MMCVIDI[3] < 0) {
                // PE
                IPC::computeRelDX_PE(mesh.vertexes[-MMCVIDI[0] - 1] - mesh.V_prev[-MMCVIDI[0] - 1],
                    mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]],
                    mesh.vertexes[MMCVIDI[2]] - mesh.V_prev[MMCVIDI[2]],
                    mesh.MMDistCoord[cI][0], relDX3D);

                Eigen::Vector2d relDX = (relDX3D * mesh.MMTanBasis[cI]).transpose();
                double relDXSqNorm = relDX.squaredNorm();
                if (relDXSqNorm > eps2) {
                    relDX /= std::sqrt(relDXSqNorm);
                }
                else {
                    double f1_div_relDXNorm;
                    IPC::f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);
                    relDX *= f1_div_relDXNorm;
                }

                Eigen::Matrix<double, 9, 1> TTTDX;
                IPC::liftRelDXTanToMesh_PE(relDX, mesh.MMTanBasis[cI], mesh.MMDistCoord[cI][0], TTTDX);
                TTTDX *= coef * mesh.Self_lambda_lastH[cI];

                grad_inc[-MMCVIDI[0] - 1] += TTTDX.template segment<3>(0);
                grad_inc[MMCVIDI[1]] += TTTDX.template segment<3>(3);
                grad_inc[MMCVIDI[2]] += TTTDX.template segment<3>(6);
            }
            else {
                // PT
                IPC::computeRelDX_PT(mesh.vertexes[-MMCVIDI[0] - 1] - mesh.V_prev[-MMCVIDI[0] - 1],
                    mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]],
                    mesh.vertexes[MMCVIDI[2]] - mesh.V_prev[MMCVIDI[2]],
                    mesh.vertexes[MMCVIDI[3]] - mesh.V_prev[MMCVIDI[3]],
                    mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], relDX3D);

                Eigen::Vector2d relDX = (relDX3D * mesh.MMTanBasis[cI]).transpose();
                double relDXSqNorm = relDX.squaredNorm();
                if (relDXSqNorm > eps2) {
                    relDX /= std::sqrt(relDXSqNorm);
                }
                else {
                    double f1_div_relDXNorm;
                    IPC::f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);
                    relDX *= f1_div_relDXNorm;
                }

                Eigen::Matrix<double, 12, 1> TTTDX;
                IPC::liftRelDXTanToMesh_PT(relDX, mesh.MMTanBasis[cI],
                    mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], TTTDX);
                TTTDX *= coef * mesh.Self_lambda_lastH[cI];

                grad_inc[-MMCVIDI[0] - 1] += TTTDX.template segment<3>(0);
                grad_inc[MMCVIDI[1]] += TTTDX.template segment<3>(3);
                grad_inc[MMCVIDI[2]] += TTTDX.template segment<3>(6);
                grad_inc[MMCVIDI[3]] += TTTDX.template segment<3>(9);
            }
        }
    }
    // ground
    int contactPairI = 0;
    for (const auto& vI : mesh.Environment_activeSet_lastH) {
        Eigen::Matrix<double, 3, 1> VDiff = (mesh.vertexes[vI] - mesh.V_prev[vI]).transpose();
        //        VDiff -= Base::velocitydt;
        Eigen::Matrix<double, 3, 1> VProj = VDiff - VDiff.dot(grd.normal) * grd.normal;
        double VProjMag2 = VProj.squaredNorm();
        if (VProjMag2 > eps2) {
            grad_inc[vI] +=
                coef * mesh.Environment_lambda_lastH[contactPairI] / std::sqrt(VProjMag2) * VProj;

        }
        else {
            grad_inc[vI] += coef * mesh.Environment_lambda_lastH[contactPairI] / eps * VProj;

        }
        ++contactPairI;
    }
}

void compute_fiction_hessian(mesh3D& mesh, const Ground& grd, BHessian& BH, double eps2, double coef) {
    double eps = std::sqrt(eps2);

    int contactPairI = 0;
    for (const auto& vI : mesh.Environment_activeSet_lastH) {
        double multiplier_vI = coef * mesh.Environment_lambda_lastH[contactPairI];
        Eigen::Matrix<double, 3, 3> H_vI;

        Eigen::Matrix<double, 3, 1> VDiff = (mesh.vertexes[vI] - mesh.V_prev[vI]).transpose();
        //        VDiff -= Base::velocitydt;
        Eigen::Matrix<double, 3, 1> VProj = VDiff - VDiff.dot(grd.normal) * grd.normal;
        double VProjMag2 = VProj.squaredNorm();
        if (VProjMag2 > eps2) {
            double VProjMag = std::sqrt(VProjMag2);

            H_vI = (VProj * (-multiplier_vI / VProjMag2 / VProjMag)) * VProj.transpose();
            H_vI += (Eigen::Matrix<double, 3, 3>::Identity() - grd.normal * grd.normal.transpose()) *
                (multiplier_vI / VProjMag);

            IglUtils::makePD<double, 3>(H_vI);
        }
        else {
            H_vI = (Eigen::Matrix<double, 3, 3>::Identity() - grd.normal * grd.normal.transpose()) *
                (multiplier_vI / eps);
            // already SPD
        }

        BH.H3x3.push_back(H_vI);
        BH.D1Index.push_back(vI);

        ++contactPairI;
    }
    // self

    std::vector<Eigen::Matrix<double, 12, 12>> IPHessian(mesh.Self_activeSet_lastH.size());
    tbb::spin_mutex countMutex2, countMutex3, countMutex4;
    tbb::parallel_for(0, (int)mesh.Self_activeSet_lastH.size(), 1, [&](int cI) {
        const auto& MMCVIDI = mesh.Self_activeSet_lastH[cI];
        if (MMCVIDI[0] >= 0) {
            // edge-edge
            Eigen::RowVector3d relDX3D;
            IPC::computeRelDX_EE(mesh.vertexes[MMCVIDI[0]] - mesh.V_prev[MMCVIDI[0]],
                mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]],
                mesh.vertexes[MMCVIDI[2]] - mesh.V_prev[MMCVIDI[2]],
                mesh.vertexes[MMCVIDI[3]] - mesh.V_prev[MMCVIDI[3]],
                mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], relDX3D);

            Eigen::Vector2d relDX = (relDX3D * mesh.MMTanBasis[cI]).transpose();
            double relDXSqNorm = relDX.squaredNorm();
            double relDXNorm = std::sqrt(relDXSqNorm);

            IPC::computeTTT_EE(mesh.MMTanBasis[cI], mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], IPHessian[cI]);
            if (relDXSqNorm > eps2) {
                IPHessian[cI] *= coef * mesh.Self_lambda_lastH[cI] / relDXNorm;

                Eigen::Matrix<double, 12, 1> TTTDX;
                IPC::liftRelDXTanToMesh_EE(relDX, mesh.MMTanBasis[cI],
                    mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], TTTDX);
                IPHessian[cI] -=
                    (TTTDX * (coef * mesh.Self_lambda_lastH[cI] / (relDXSqNorm * relDXNorm))) * TTTDX.transpose();

                IglUtils::makePD<double, 12>(IPHessian[cI]);
            }
            else {
                double f1_div_relDXNorm;
                IPC::f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);

                IPHessian[cI] *= coef * mesh.Self_lambda_lastH[cI] * f1_div_relDXNorm;

                double f2;
                IPC::f2_SF(relDXSqNorm, eps, f2);
                if (f2 != f1_div_relDXNorm && relDXSqNorm) {
                    // higher order clamping and not exactly static
                    Eigen::Matrix<double, 12, 1> TTTDX;
                    IPC::liftRelDXTanToMesh_EE(relDX, mesh.MMTanBasis[cI],
                        mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], TTTDX);
                    IPHessian[cI] +=
                        (TTTDX * (coef * mesh.Self_lambda_lastH[cI] * (f2 - f1_div_relDXNorm) / relDXSqNorm)) *
                        TTTDX.transpose();

                    IglUtils::makePD<double, 12>(IPHessian[cI]);
                }
            }
            countMutex4.lock();
            BH.H12x12.emplace_back(IPHessian[cI]);
            BH.D4Index.emplace_back(MMCVIDI[0], MMCVIDI[1], MMCVIDI[2], MMCVIDI[3]);
            countMutex4.unlock();
        }
        else {
            // point-triangle and degenerate edge-edge
            assert(MMCVIDI[1] >= 0);
            int v0I = -MMCVIDI[0] - 1;
            if (MMCVIDI[2] < 0) {
                // PP
                Eigen::RowVector3d relDX3D;
                IPC::computeRelDX_PP(mesh.vertexes[v0I] - mesh.V_prev[v0I],
                    mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]], relDX3D);

                Eigen::Vector2d relDX = (relDX3D * mesh.MMTanBasis[cI]).transpose();
                double relDXSqNorm = relDX.squaredNorm();
                double relDXNorm = std::sqrt(relDXSqNorm);

                Eigen::Matrix<double, 6, 6> HessianBlock;
                IPC::computeTTT_PP(mesh.MMTanBasis[cI], HessianBlock);
                if (relDXSqNorm > eps2) {
                    HessianBlock *= coef * mesh.Self_lambda_lastH[cI] / relDXNorm;

                    Eigen::Matrix<double, 6, 1> TTTDX;
                    IPC::liftRelDXTanToMesh_PP(relDX, mesh.MMTanBasis[cI], TTTDX);
                    HessianBlock -=
                        (TTTDX * (coef * mesh.Self_lambda_lastH[cI] / (relDXSqNorm * relDXNorm))) *
                        TTTDX.transpose();

                    IglUtils::makePD<double, 6>(HessianBlock);
                }
                else {
                    double f1_div_relDXNorm;
                    IPC::f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);

                    HessianBlock *= coef * mesh.Self_lambda_lastH[cI] * f1_div_relDXNorm;

                    double f2;
                    IPC::f2_SF(relDXSqNorm, eps, f2);
                    if (f2 != f1_div_relDXNorm && relDXSqNorm) {
                        // higher order clamping and not exactly static
                        Eigen::Matrix<double, 6, 1> TTTDX;
                        IPC::liftRelDXTanToMesh_PP(relDX, mesh.MMTanBasis[cI], TTTDX);
                        HessianBlock +=
                            (TTTDX * (coef * mesh.Self_lambda_lastH[cI] * (f2 - f1_div_relDXNorm) / relDXSqNorm)) *
                            TTTDX.transpose();

                        IglUtils::makePD<double, 6>(HessianBlock);
                    }
                }
                countMutex2.lock();
                BH.H6x6.emplace_back(HessianBlock);
                BH.D2Index.emplace_back(v0I, MMCVIDI[1]);
                countMutex2.unlock();
            }
            else if (MMCVIDI[3] < 0) {
                // PE
                Eigen::RowVector3d relDX3D;
                IPC::computeRelDX_PE(mesh.vertexes[v0I] - mesh.V_prev[v0I],
                    mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]],
                    mesh.vertexes[MMCVIDI[2]] - mesh.V_prev[MMCVIDI[2]],
                    mesh.MMDistCoord[cI][0], relDX3D);

                Eigen::Vector2d relDX = (relDX3D * mesh.MMTanBasis[cI]).transpose();
                double relDXSqNorm = relDX.squaredNorm();
                double relDXNorm = std::sqrt(relDXSqNorm);

                Eigen::Matrix<double, 9, 9> HessianBlock;
                IPC::computeTTT_PE(mesh.MMTanBasis[cI], mesh.MMDistCoord[cI][0], HessianBlock);
                if (relDXSqNorm > eps2) {
                    HessianBlock *= coef * mesh.Self_lambda_lastH[cI] / relDXNorm;

                    Eigen::Matrix<double, 9, 1> TTTDX;
                    IPC::liftRelDXTanToMesh_PE(relDX, mesh.MMTanBasis[cI], mesh.MMDistCoord[cI][0], TTTDX);
                    HessianBlock -=
                        (TTTDX * (coef * mesh.Self_lambda_lastH[cI] / (relDXSqNorm * relDXNorm))) *
                        TTTDX.transpose();

                    IglUtils::makePD<double, 9>(HessianBlock);
                }
                else {
                    double f1_div_relDXNorm;
                    IPC::f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);

                    HessianBlock *= coef * mesh.Self_lambda_lastH[cI] * f1_div_relDXNorm;

                    double f2;
                    IPC::f2_SF(relDXSqNorm, eps, f2);
                    if (f2 != f1_div_relDXNorm && relDXSqNorm) {
                        // higher order clamping and not exactly static
                        Eigen::Matrix<double, 9, 1> TTTDX;
                        IPC::liftRelDXTanToMesh_PE(relDX, mesh.MMTanBasis[cI], mesh.MMDistCoord[cI][0], TTTDX);
                        HessianBlock +=
                            (TTTDX * (coef * mesh.Self_lambda_lastH[cI] * (f2 - f1_div_relDXNorm) / relDXSqNorm)) *
                            TTTDX.transpose();

                        IglUtils::makePD<double, 9>(HessianBlock);
                    }
                }
                countMutex3.lock();
                BH.H9x9.emplace_back(HessianBlock);
                BH.D3Index.emplace_back(v0I, MMCVIDI[1], MMCVIDI[2]);
                countMutex3.unlock();
            }
            else {
                // PT
                Eigen::RowVector3d relDX3D;
                IPC::computeRelDX_PT(mesh.vertexes[v0I] - mesh.V_prev[v0I],
                    mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]],
                    mesh.vertexes[MMCVIDI[2]] - mesh.V_prev[MMCVIDI[2]],
                    mesh.vertexes[MMCVIDI[3]] - mesh.V_prev[MMCVIDI[3]],
                    mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], relDX3D);

                Eigen::Vector2d relDX = (relDX3D * mesh.MMTanBasis[cI]).transpose();
                double relDXSqNorm = relDX.squaredNorm();
                double relDXNorm = std::sqrt(relDXSqNorm);

                IPC::computeTTT_PT(mesh.MMTanBasis[cI], mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1],
                    IPHessian[cI]);
                if (relDXSqNorm > eps2) {
                    IPHessian[cI] *= coef * mesh.Self_lambda_lastH[cI] / relDXNorm;
                    Eigen::Matrix<double, 12, 1> TTTDX;
                    IPC::liftRelDXTanToMesh_PT(relDX, mesh.MMTanBasis[cI],
                        mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], TTTDX);
                    IPHessian[cI] -=
                        (TTTDX * (coef * mesh.Self_lambda_lastH[cI] / (relDXSqNorm * relDXNorm))) *
                        TTTDX.transpose();

                    IglUtils::makePD<double, 12>(IPHessian[cI]);
                }
                else {
                    double f1_div_relDXNorm;
                    IPC::f1_SF_div_relDXNorm(relDXSqNorm, eps, f1_div_relDXNorm);

                    IPHessian[cI] *= coef * mesh.Self_lambda_lastH[cI] * f1_div_relDXNorm;

                    double f2;
                    IPC::f2_SF(relDXSqNorm, eps, f2);
                    if (f2 != f1_div_relDXNorm && relDXSqNorm) {
                        // higher order clamping and not exactly static
                        Eigen::Matrix<double, 12, 1> TTTDX;
                        IPC::liftRelDXTanToMesh_PT(relDX, mesh.MMTanBasis[cI],
                            mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], TTTDX);
                        IPHessian[cI] +=
                            (TTTDX * (coef * mesh.Self_lambda_lastH[cI] * (f2 - f1_div_relDXNorm) / relDXSqNorm)) *
                            TTTDX.transpose();

                        IglUtils::makePD<double, 12>(IPHessian[cI]);
                    }
                }
                countMutex4.lock();
                BH.H12x12.emplace_back(IPHessian[cI]);
                BH.D4Index.emplace_back(v0I, MMCVIDI[1], MMCVIDI[2], MMCVIDI[3]);
                countMutex4.unlock();
            }
        }
        }
    );


}

void computeXTilta(mesh3D& mesh) {
    mesh.xTilta.resize(mesh.vertexes.size());
	double gravity = SimulationParameters::gravityMagnitude;
    if(!mesh.apply_gravity) {
        gravity = 0;
	}
    Vector3d gravityDtSq = Vector3d(0, gravity, 0) * mesh.IPC_dt * mesh.IPC_dt;

    tbb::parallel_for(0, (int)mesh.vertexes.size(), 1, [&](int vI) {
        mesh.xTilta[vI] = (mesh.V_prev[vI] + (mesh.velocities[vI] * mesh.IPC_dt + gravityDtSq));
        }

    );

}

void computeEGradient(const mesh3D& mesh, vector<Vector3d>& gradient) {
    //Mass part
    Matrix3d massM;
    massM.setIdentity();

    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int vI) {
        gradient[vI] += (mesh.masses[vI] * massM * (mesh.vertexes[vI] - mesh.xTilta[vI]));
        }

    );

    //internal part
    //double tolerance = 1e-6;
    double fsum = 0;

    vector<tbb::spin_mutex> countMutex(mesh.vertexNum);
    tbb::parallel_for(0, mesh.tetrahedraNum, 1, [&](int ii) {
        const Eigen::Matrix<double, 9, 12>& PFPX = mesh.tetPFPX[ii];
        MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_tetrahedra_inverse[ii];

        Matrix3d PEPF = computePEPF_StableNHK3D_2_double(F, mesh.lengthRate, mesh.volumeRate);

        Eigen::Map<const Eigen::Matrix<double, 9, 1>> pepf(PEPF.data());
        Eigen::Matrix<double, 12, 1> f = mesh.volum[ii] * PFPX.transpose() * pepf;


        for (int i = 0; i < 12; i++) {

            countMutex[mesh.tetrahedras[ii][i / 3]].lock();
            gradient[mesh.tetrahedras[ii][i / 3]][i % 3] += mesh.IPC_dt * mesh.IPC_dt * f(i, 0);
            countMutex[mesh.tetrahedras[ii][i / 3]].unlock();

        }
        }

    );
}

static void addInertiaAndDragGradient(mesh3D& mesh, vector<Vector3d>& gradient) {
    //calculate inertial gradient
    Matrix3d massM;
    massM.setIdentity();
    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int vI) {
        gradient[vI] = (mesh.masses[vI] * massM * (mesh.vertexes[vI] - mesh.xTilta[vI]));
        }

    );
    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int vI) {
        gradient[vI] += (mesh.drag_coeff * mesh.masses[vI] * massM * (mesh.vertexes[vI] - mesh.V_prev[vI]));
        }

    );
}

static void addTetrahedralElasticityGradientHessian(mesh3D& mesh, vector<Vector3d>& gradient, BHessian& BH) {
    vector<tbb::spin_mutex> countMutex(mesh.vertexNum);
    tbb::parallel_for(0, mesh.tetrahedraNum, 1, [&](int ii) {
        const Eigen::Matrix<double, 9, 12>& PFPX = mesh.tetPFPX[ii];
        MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_tetrahedra_inverse[ii];
        //cout << mesh.lengthRate << "   " << mesh.volumeRate << endl;
        Matrix3d PEPF = computePEPF_StableNHK3D_2_double(F, mesh.lengthRate, mesh.volumeRate);

        Eigen::Map<const Eigen::Matrix<double, 9, 1>> pepf(PEPF.data());
        Eigen::Matrix<double, 12, 1> f = mesh.volum[ii] * PFPX.transpose() * pepf;

        for (int i = 0; i < 4; i++) {
            countMutex[mesh.tetrahedras[ii][i]].lock();
            gradient[mesh.tetrahedras[ii][i]] += mesh.IPC_dt * mesh.IPC_dt * f.template segment<3>(i * 3);//f(i, 0);
            countMutex[mesh.tetrahedras[ii][i]].unlock();

        }

        //std::cout << F << std::endl;
        MatrixXd Hq = project_StabbleNHK_2_H_3D(F, mesh.lengthRate, mesh.volumeRate);

        MatrixXd HE = mesh.volum[ii] * mesh.IPC_dt * mesh.IPC_dt * PFPX.transpose() * Hq * PFPX;
        BH.H12x12[ii] = (HE);

        }

    );
    BH.D4Index = mesh.tetrahedras;
}

static void addTriangleShellGradientHessian(mesh3D& mesh, vector<Vector3d>& gradient, BHessian& BH) {
    vector<tbb::spin_mutex> countMutex(mesh.vertexNum);
    BH.H9x9.resize(mesh.triangleNum);
    tbb::parallel_for(0, mesh.triangleNum, 1, [&](int ii) {
        const Eigen::Matrix<double, 6, 9>& PFPX = mesh.triPFPX[ii];
        Matrix<double, 3, 2> F =
            calculateDs32D_double(mesh.vertexes, mesh.triangles[ii]) * mesh.DM_triangle_inverse[ii];
        Vector2d anisotropic_a = Vector2d(1, 0), anisotropic_b = Vector2d(0, 1);
        Matrix<double, 3, 2> PEPF = computePEPF_baraffwitkin_double(F, anisotropic_a, anisotropic_b, mesh.stretchStiffness,
            mesh.shearStiffness, mesh.strainRate);

        Eigen::Map<const Eigen::Matrix<double, 6, 1>> pepf(PEPF.data());
        Eigen::Matrix<double, 9, 1> f = mesh.IPC_dt * mesh.IPC_dt * mesh.areas[ii] * PFPX.transpose() * pepf;

        for (int i = 0; i < 3; i++) {
            countMutex[mesh.triangles[ii][i]].lock();
            gradient[mesh.triangles[ii][i]] += f.template segment<3>(i * 3);//f(i, 0);
            countMutex[mesh.triangles[ii][i]].unlock();
        }

        //std::cout << F << std::endl;
        MatrixXd Hq = project_baraffwitkint_H_3D(F, anisotropic_a, anisotropic_b, mesh.stretchStiffness,
            mesh.shearStiffness, mesh.strainRate);

        MatrixXd HE = mesh.areas[ii] * mesh.IPC_dt * mesh.IPC_dt * PFPX.transpose() * Hq * PFPX;
        BH.H9x9[ii] = (HE);

        }

    );
    BH.D3Index = mesh.triangles;
}

static void addBendingGradientHessian(mesh3D& mesh, vector<Vector3d>& gradient, BHessian& BH) {
    vector<tbb::spin_mutex> countMutex(mesh.vertexNum);

#if defined USE_QUADRATIC_BENDING
    size_t bendingInfoSize = mesh.quadBendingInfo.size();
    std::vector<Vector3d>& verts = mesh.vertexes;

    vector<Vector4i> bendIndexes(bendingInfoSize);

    // gradient
    tbb::parallel_for(
        0, static_cast<int>(bendingInfoSize), 1, [&](int i)
        {
            const QuadBendingInfo& info = mesh.quadBendingInfo[i];

            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    for (int d = 0; d < 3; ++d)
                    {
                        countMutex[info.verts[i]].lock();
                        double incr = mesh.bendingStiffness * mesh.IPC_dt * mesh.IPC_dt * info.Q(i, j) * verts[info.verts[j]][d];
                        gradient[info.verts[i]][d] += incr;
                        //bendingGradient[info.verts[i]][d] += incr;
                        countMutex[info.verts[i]].unlock();
                    }
        }
    );

    // hessian
    size_t H12Offset = mesh.tetrahedraNum;

    tbb::parallel_for(
        0, static_cast<int>(bendingInfoSize), 1, [&](int i)
        {
            const QuadBendingInfo& info = mesh.quadBendingInfo[i];
            Matrix<double, 12, 12> H = Matrix<double, 12, 12>::Zero();
            for (int i = 0; i != 4; ++i)
                for (int j = 0; j != 4; ++j)
                    for (int d = 0; d != 3; ++d)
                        H(3 * i + d, 3 * j + d) = mesh.bendingStiffness * mesh.IPC_dt * mesh.IPC_dt * info.Q(i, j);

            Eigen::EigenSolver<Matrix<double, -1, -1>> es(H);

            Eigen::VectorXd eigenvalues = es.eigenvalues().real();

            for (int i = 0; i != eigenvalues.size(); ++i)
                if (eigenvalues[i] < 0)
                    eigenvalues[i] = 0;

            Eigen::MatrixXd V = es.eigenvectors().real();
            Eigen::MatrixXd Vt = V.transpose();
            Eigen::MatrixXd lambda = eigenvalues.asDiagonal();

            H = V * lambda * Vt;

            Vector4i localIndex = Vector4i(info.verts[0], info.verts[1], info.verts[2], info.verts[3]);

            BH.H12x12[H12Offset + i] = H;
            bendIndexes[i] = localIndex;
        }
    );

    BH.D4Index.insert(BH.D4Index.end(), bendIndexes.begin(), bendIndexes.end());
#else
    int offset = mesh.tetrahedraNum;
    vector<Vector4i> bendIndexes(mesh.tri_edges.size());
    //bending
    //for (int idx = 0;idx < mesh.tri_edges.size();idx++)
    tbb::parallel_for(0, (int)mesh.tri_edges.size(), 1, [&](int idx)
        {
            const auto& edge = mesh.tri_edges[idx];
            const auto& adj_verts = mesh.tri_edges_adj_points[idx];
            if (adj_verts.y() == -1) {
                Matrix<double, 12, 12> zeroM;
                zeroM.setZero();

                BH.H12x12[offset + idx] = zeroM;
                bendIndexes[idx] = Vector4i(0, 1, 2, 3);

            }
            else {
                auto x0 = mesh.vertexes[edge.x()];
                auto x1 = mesh.vertexes[edge.y()];
                auto x2 = mesh.vertexes[adj_verts.x()];
                auto x3 = mesh.vertexes[adj_verts.y()];

                Matrix<double, 1, 12> grad_transpose;
                Matrix<double, 12, 12> H;

                double t = edgeTheta(x0, x1, x2, x3, &grad_transpose, &H);
                //            cout << "t: " << t << endl;

                auto rest_x0 = mesh.v_rest[edge.x()];
                auto rest_x1 = mesh.v_rest[edge.y()];
                auto rest_x2 = mesh.v_rest[adj_verts.x()];
                auto rest_x3 = mesh.v_rest[adj_verts.y()];
                double length = (rest_x0 - rest_x1).norm();

                double rest_t = edgeTheta(rest_x0, rest_x1, rest_x2, rest_x3, nullptr, nullptr);

                H = 2 * ((t - rest_t) * H + grad_transpose.transpose() * grad_transpose);
                grad_transpose = 2 * (t - rest_t) * grad_transpose;
                IglUtils::makePD<double, 12>(H);

                VectorXd f = mesh.bendingStiffness * mesh.IPC_dt * mesh.IPC_dt * length * grad_transpose;

                for (int i = 0; i < 2; i++) {
                    countMutex[edge[i]].lock();
                    gradient[edge[i]] += f.template segment<3>(i * 3);//f(i, 0);
                    countMutex[edge[i]].unlock();

                    countMutex[adj_verts[i]].lock();
                    gradient[adj_verts[i]] += f.template segment<3>((i + 2) * 3);
                    countMutex[adj_verts[i]].unlock();

                }

                BH.H12x12[offset + idx] = (mesh.bendingStiffness * mesh.IPC_dt * mesh.IPC_dt * length * H);
                bendIndexes[idx] = (Vector4i(edge.x(), edge.y(), adj_verts.x(), adj_verts.y()));
            }

            //bendIndexes[idx] = Vector4i(edge.x(), edge.y(), adj_verts.x(), adj_verts.y());
        }
    );

    BH.D4Index.insert(BH.D4Index.end(), bendIndexes.begin(), bendIndexes.end());
#endif
}

static void addGroundBarrierGradientHessian(mesh3D& mesh, const Ground& grd, vector<Vector3d>& gradient, BHessian& BH) {
    //calculate barrier gradient and hessian for the ground plane collision
    VectorXd constraintVals;
    int offset = 0;
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
}

static int addSelfContactGradientHessian(mesh3D& mesh, vector<Vector3d>& gradient, BHessian& BH) {
    int collisionNum = 0;

    VectorXd constraintVals;
    Evaluate_SelfPTConstraintVals(mesh, constraintVals, 0);

    tbb::parallel_for(0, (int)constraintVals.size(), 1, [&](int cI)
        //for (int cI = offset; cI < constraintVals.size(); ++cI)
        {
            compute_g_b(constraintVals[cI], mesh.Hhat, constraintVals[cI]);
        }
    );
#ifndef NEWB
    collisionNum += compute_g_dpt(mesh, mesh.Self_ActiveSet, constraintVals, gradient, 0, mesh.Kappa);
#else
    compute_g_dpt_new(mesh, mesh.Self_ActiveSet, constraintVals, gradient, 0, mesh.Kappa, mesh.Hhat);
#endif

    compute_g_dee(mesh, gradient, mesh.Hhat, mesh.Kappa);

#ifndef NEWB
    compute_H_dpt(mesh, BH, mesh.Hhat, mesh.Kappa);
#else
    compute_H_dpt_new(mesh, BH, mesh.Hhat, mesh.Kappa, mesh.Hhat);
#endif
    compute_H_dee(mesh, BH, mesh.Hhat, mesh.Kappa);

    return collisionNum;
}

static void addFrictionGradientHessian(mesh3D& mesh, const Ground& grd, vector<Vector3d>& gradient, BHessian& BH) {
#ifdef USE_FRICTION
    compute_fiction_gradient(mesh, grd, gradient, mesh.Fhat * mesh.IPC_dt * mesh.IPC_dt, mesh.friction);
    compute_fiction_hessian(mesh, grd, BH, mesh.Fhat * mesh.IPC_dt * mesh.IPC_dt, mesh.friction);
#endif
}

int computeGradientAndHessian(mesh3D& mesh, vector<Vector3d>& gradient, BHessian& BH, const Ground& grd) {
    addInertiaAndDragGradient(mesh, gradient);

#if defined USE_QUADRATIC_BENDING
    BH.H12x12.resize(mesh.tetrahedraNum + mesh.quadBendingInfo.size());
#else
    BH.H12x12.resize(mesh.tetrahedraNum + mesh.tri_edges.size());
#endif

    addTetrahedralElasticityGradientHessian(mesh, gradient, BH);
    addTriangleShellGradientHessian(mesh, gradient, BH);
    addBendingGradientHessian(mesh, gradient, BH);
    addGroundBarrierGradientHessian(mesh, grd, gradient, BH);
    int collisionNum = addSelfContactGradientHessian(mesh, gradient, BH);
    addFrictionGradientHessian(mesh, grd, gradient, BH);

    return collisionNum;
}


void __Inverse2(const Eigen::Matrix3d& input, Eigen::Matrix3d& result) {
    double eps = 1e-15;
    const int dim = 3;
    double mat[dim * 2][dim * 2];
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < 2 * dim; j++) {
            if (j < dim) {
                mat[i][j] = input(i, j);
            }
            else {
                mat[i][j] = j - dim == i ? 1 : 0;
            }
        }
    }

    for (int i = 0; i < dim; i++) {
        if (abs(mat[i][i]) < eps) {
            int j;
            for (j = i + 1; j < dim; j++) {
                if (abs(mat[j][i]) > eps) break;
            }
            if (j == dim) return;
            for (int r = i; r < 2 * dim; r++) {
                mat[i][r] += mat[j][r];
            }
        }
        double ep = mat[i][i];
        for (int r = i; r < 2 * dim; r++) {
            mat[i][r] /= ep;
        }

        for (int j = i + 1; j < dim; j++) {
            double e = -1 * (mat[j][i] / mat[i][i]);
            for (int r = i; r < 2 * dim; r++) {
                mat[j][r] += e * mat[i][r];
            }
        }
    }

    for (int i = dim - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            double e = -1 * (mat[j][i] / mat[i][i]);
            for (int r = i; r < 2 * dim; r++) {
                mat[j][r] += e * mat[i][r];
            }
        }
    }


    for (int i = 0; i < dim; i++) {
        for (int r = dim; r < 2 * dim; r++) {
            result(i, r - dim) = mat[i][r];
        }
    }
}

void PCG_Precondition(const mesh3D& mesh, const BHessian& BH, const vector<Vector3d>& gradient, vector<Matrix3d>& P) {


    vector<tbb::spin_mutex> countMutexVerts(mesh.vertexNum);//, countMutex9, countMutex6;//, countMutex3;
    tbb::parallel_for(0, (int)mesh.vertexNum, 1, [&](int i) {
        P[i].setZero();
        P[i](0, 0) = mesh.masses[i];
        P[i](1, 1) = mesh.masses[i];
        P[i](2, 2) = mesh.masses[i];
        });


    //for (int i = 0;i < BH.D4Index.size();i++) 
    tbb::parallel_for(0, (int)BH.D4Index.size(), 1, [&](int i) {
        for (int j = 0; j < 12; j++) {
            //P[BH.D4Index[i][j / 3]][j % 3] += BH.H12x12[i](j, j);
            countMutexVerts[BH.D4Index[i][j / 3]].lock();
            P[BH.D4Index[i][j / 3]](0, j % 3) += BH.H12x12[i]((j / 3) * 3, j);
            P[BH.D4Index[i][j / 3]](1, j % 3) += BH.H12x12[i]((j / 3) * 3 + 1, j);
            P[BH.D4Index[i][j / 3]](2, j % 3) += BH.H12x12[i]((j / 3) * 3 + 2, j);
            countMutexVerts[BH.D4Index[i][j / 3]].unlock();
        }
        });
    //for (int i = 0;i < BH.D3Index.size();i++) 
    tbb::parallel_for(0, (int)BH.D3Index.size(), 1, [&](int i) {
        for (int j = 0; j < 9; j++) {
            //P[BH.D3Index[i][j / 3]][j % 3] += BH.H9x9[i](j, j);
            countMutexVerts[BH.D3Index[i][j / 3]].lock();
            P[BH.D3Index[i][j / 3]](0, j % 3) += BH.H9x9[i]((j / 3) * 3, j);
            P[BH.D3Index[i][j / 3]](1, j % 3) += BH.H9x9[i]((j / 3) * 3 + 1, j);
            P[BH.D3Index[i][j / 3]](2, j % 3) += BH.H9x9[i]((j / 3) * 3 + 2, j);
            countMutexVerts[BH.D3Index[i][j / 3]].unlock();
        }
        });
    //for (int i = 0;i < BH.D2Index.size();i++) 
    tbb::parallel_for(0, (int)BH.D2Index.size(), 1, [&](int i) {
        for (int j = 0; j < 6; j++) {
            //P[BH.D2Index[i][j / 3]][j % 3] += BH.H6x6[i](j, j);
            countMutexVerts[BH.D2Index[i][j / 3]].lock();
            P[BH.D2Index[i][j / 3]](0, j % 3) += BH.H6x6[i]((j / 3) * 3, j);
            P[BH.D2Index[i][j / 3]](1, j % 3) += BH.H6x6[i]((j / 3) * 3 + 1, j);
            P[BH.D2Index[i][j / 3]](2, j % 3) += BH.H6x6[i]((j / 3) * 3 + 2, j);
            countMutexVerts[BH.D2Index[i][j / 3]].unlock();
        }
        });
    //for (int i = 0;i < BH.D1Index.size();i++) 
    tbb::parallel_for(0, (int)BH.D1Index.size(), 1, [&](int i) {
        for (int j = 0; j < 3; j++) {
            //P[BH.D1Index[i]][j] += BH.H3x3[i](j, j);
            countMutexVerts[BH.D1Index[i]].lock();
            P[BH.D1Index[i]](0, j) += BH.H3x3[i](0, j);
            P[BH.D1Index[i]](1, j) += BH.H3x3[i](1, j);
            P[BH.D1Index[i]](2, j) += BH.H3x3[i](2, j);
            countMutexVerts[BH.D1Index[i]].unlock();
        }
        });
    //for (int i = 0;i < mesh.vertexNum;i++) 
    tbb::parallel_for(0, (int)mesh.vertexNum, 1, [&](int i) {
        Eigen::Matrix3d DMInverse;
        __Inverse2(P[i], DMInverse);

        P[i] = DMInverse;
        });
}

vector<Vector3d> cholmod_solver_EPF(BHessian& BH, std::vector<Eigen::Vector3d> gradient, const mesh3D& mesh, vector<Eigen::VectorXd>& KK) {

    vector<Vector3d> direction(mesh.vertexNum, Vector3d(0, 0, 0));
    const int vectorSize = 3 * mesh.vertexNum;
    CholmodSolver solver;
    std::vector<Triplet<double>> triplets = BH.toTriplets(mesh.boundaryTypes);
    int offset = triplets.size();
    triplets.resize(offset + vectorSize);
    tbb::parallel_for(0, vectorSize, 1, [&](int i) {
        triplets[offset + i] = Triplet<double>(i, i, mesh.masses[i / 3]);
        }
    );
    Eigen::SparseMatrix<double> sparseFEM(vectorSize, vectorSize);
    sparseFEM.setFromTriplets(triplets.begin(), triplets.end());

    solver.preFactorize(sparseFEM);


    auto vec_gradient = Eigen::VectorXd(vectorSize);

    tbb::parallel_for(0, (int)(gradient.size()), 1, [&](int i) {
        if (mesh.boundaryTypes[i] == 0)
            vec_gradient.block<3, 1>(3 * i, 0) = gradient[i];
        else
            vec_gradient.block<3, 1>(3 * i, 0) = Vector3d(0, 0, 0);
        }
    );
    auto result0 = Eigen::VectorXd(gradient.size() * 3);

    solver.solve_with_preFactorize(vec_gradient, result0);
    tbb::parallel_for(0, (int)(gradient.size()), 1, [&](int i) {
        direction[i] = result0.block<3, 1>(3 * i, 0);
        }
    );

    for (int i = 0; i < mesh.vertexNum; i++) {
        for (int j = 0; j < 3; j++) {
            vec_gradient.setZero();
            vec_gradient(i * 3 + j) = mesh.masses[i];
            auto result = Eigen::VectorXd(vectorSize);
            solver.solve_with_preFactorize(vec_gradient, result);
            KK.push_back(result);
        }
    }

    return direction;
}


vector<Vector3d> Eigen_CG_solver(BHessian& BH,
    std::vector<Eigen::Vector3d> gradient,
    const mesh3D& mesh)
{

    vector<Vector3d> direction(mesh.vertexNum, Vector3d(0, 0, 0));


    //printf("start triplets\n");
    std::vector<Triplet<double>> triplets =
        BH.toTriplets(mesh.boundaryTypes);
    //printf("finish triplets\n");
    int offset = triplets.size();
    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //if (!mesh.is_quasi_static)
    {
        triplets.resize(offset + 3 * mesh.vertexNum);

        tbb::parallel_for(0,
            3 * mesh.vertexNum,
            1,
            [&](int i) {
                triplets[offset + i] =
                    Triplet<double>(i, i, mesh.masses[i / 3]);
            });
    }

    Eigen::SparseMatrix<double> sparseFEM(3 * mesh.vertexNum, 3 * mesh.vertexNum);
    sparseFEM.setFromTriplets(triplets.begin(), triplets.end());
    auto vec_gradient = Eigen::VectorXd(gradient.size() * 3);

    tbb::parallel_for(0,
        (int)(gradient.size()),
        1,
        [&](int i)
        {
            if (mesh.boundaryTypes[i] == 0)
                vec_gradient.block<3, 1>(3 * i, 0) = gradient[i];
            else
                vec_gradient.block<3, 1>(3 * i, 0) = Vector3d(0, 0, 0);
        });
    auto result = Eigen::VectorXd(gradient.size() * 3);

    //Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower, IncompleteCholesky> cg;
    Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower, Eigen::IncompleteCholesky<double, Lower>> cg;
    cg.setTolerance(1e-2);
    cg.compute(sparseFEM);
    result = cg.solve(vec_gradient);

    tbb::parallel_for(0,
        (int)(gradient.size()),
        1,
        [&](int i) { direction[i] = result.block<3, 1>(3 * i, 0); });
    //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //std::cout << "Time for computing cg only = "<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;


    return direction;
}


vector<Vector3d> cholmod_solver(BHessian& BH, std::vector<Eigen::Vector3d> gradient, const mesh3D& mesh) {

    vector<Vector3d> direction(mesh.vertexNum, Vector3d(0, 0, 0));
    CholmodSolver solver;
    std::vector<Triplet<double>> triplets = BH.toTriplets(mesh.boundaryTypes);
    int offset = triplets.size();
    triplets.resize(offset + 3 * mesh.vertexNum);

    tbb::parallel_for(0, 3 * mesh.vertexNum, 1, [&](int i) {
        triplets[offset + i] = Triplet<double>(i, i, mesh.masses[i / 3] * (1 + mesh.drag_coeff));
        }
    );
    Eigen::SparseMatrix<double> sparseFEM(3 * mesh.vertexNum, 3 * mesh.vertexNum);
    sparseFEM.setFromTriplets(triplets.begin(), triplets.end());
    auto vec_gradient = Eigen::VectorXd(gradient.size() * 3);

    tbb::parallel_for(0, (int)(gradient.size()), 1, [&](int i) {
        if (mesh.boundaryTypes[i] == 0)
            vec_gradient.block<3, 1>(3 * i, 0) = gradient[i];
        else
            vec_gradient.block<3, 1>(3 * i, 0) = Vector3d(0, 0, 0);
        }
    );
    auto result = Eigen::VectorXd(gradient.size() * 3);

    solver.set_pattern(sparseFEM);
    solver.solve(vec_gradient, result);
    tbb::parallel_for(0, (int)(gradient.size()), 1, [&](int i) {
        direction[i] = result.block<3, 1>(3 * i, 0);
        }
    );

    return direction;
}

vector<Vector3d> PCG_Solver(const mesh3D& mesh, const BHessian& BH, const vector<Vector3d>& gradient) {
    //std::vector<Vector3d> P(mesh.vertexNum, Vector3d(0, 0, 0));
    std::vector<Matrix3d> P(mesh.vertexNum);
    PCG_Precondition(mesh, BH, gradient, P);

    vector<Vector3d> tempDeltaX(mesh.vertexNum, Vector3d(0, 0, 0));

    double deltaN = 0;
    double localOptimal = 1e32;//std::numeric_limits<double>::max();
    bool getLocalOpt = false;
    vector<Vector3d> dX(mesh.vertexNum, Vector3d(0, 0, 0));

    double delta0 = 0;
    double deltaO = 0;

    vector<Vector3d> r(mesh.vertexNum, Vector3d(0, 0, 0));
    vector<Vector3d> c(mesh.vertexNum, Vector3d(0, 0, 0));

    vector<tbb::spin_mutex> countMutex(mesh.vertexNum);
    delta0 = parallel_reduce(
        tbb::blocked_range<int>(0, mesh.vertexNum), 0.0, [&](const tbb::blocked_range<int>& r, double temp_delta) {
            for (int i = r.begin(); i != r.end(); i++) {

                Vector3d filter_b = mesh.Constraints[i] * gradient[i];
                temp_delta += filter_b.dot((P[i] * filter_b));
            }
            return temp_delta;
        },
        [&](double left, double right) {
            return left + right;
        }
    );

    deltaN = parallel_reduce(
        tbb::blocked_range<int>(0, mesh.vertexNum), 0.0,
        [&](const tbb::blocked_range<int>& rg, double temp_deltaN) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                r[i] = gradient[i];
                r[i] = mesh.Constraints[i] * r[i];
                c[i] = P[i] * r[i];
                c[i] = mesh.Constraints[i] * c[i];

                temp_deltaN += c[i].dot(r[i]);

            }
            return temp_deltaN;
        },
        [&](double left, double right) {
            return left + right;
        }
    );

    double errorRate = 1e-4;
    //std::cout << "cpu  delta0:   " << delta0 << "      deltaN:   " << deltaN << endl;
    //system("pause");
    //PCG main loop
    int cgCounts = 0;
    while (cgCounts < 30000 && deltaN > errorRate * delta0) {
        cgCounts++;
        //std::cout << "delta0:   " << delta0 << "      deltaN:   " << deltaN  << endl;
        vector<Vector3d> q(mesh.vertexNum, Vector3d(0, 0, 0));

        tbb::parallel_for(0, (int)BH.D4Index.size(), 1, [&](int ii) {
            MatrixXd H = BH.H12x12[ii];
            VectorXd tempC(12);
            for (int i = 0; i < 12; i++) {
                tempC(i) = c[BH.D4Index[ii][i / 3]][i % 3];
            }
            VectorXd tempQ = H * tempC;
            for (int i = 0; i < 12; i++) {
                countMutex[BH.D4Index[ii][i / 3]].lock();
                q[BH.D4Index[ii][i / 3]][i % 3] += tempQ(i);
                countMutex[BH.D4Index[ii][i / 3]].unlock();
            }
            }
        );

        tbb::parallel_for(0, (int)BH.D3Index.size(), 1, [&](int ii) {
            MatrixXd H = BH.H9x9[ii];

            VectorXd tempC(9);

            for (int i = 0; i < 9; i++) {
                tempC(i) = c[BH.D3Index[ii][i / 3]][i % 3];
            }
            VectorXd tempQ = H * tempC;

            for (int i = 0; i < 9; i++) {
                countMutex[BH.D3Index[ii][i / 3]].lock();
                q[BH.D3Index[ii][i / 3]][i % 3] += tempQ(i);
                countMutex[BH.D3Index[ii][i / 3]].unlock();
            }
            }
        );


        tbb::parallel_for(0, (int)BH.D2Index.size(), 1, [&](int ii) {
            MatrixXd H = BH.H6x6[ii];

            VectorXd tempC(6);

            for (int i = 0; i < 6; i++) {
                tempC(i) = c[BH.D2Index[ii][i / 3]][i % 3];
            }
            VectorXd tempQ = H * tempC;

            for (int i = 0; i < 6; i++) {
                countMutex[BH.D2Index[ii][i / 3]].lock();
                q[BH.D2Index[ii][i / 3]][i % 3] += tempQ(i);
                countMutex[BH.D2Index[ii][i / 3]].unlock();
            }
            }
        );

        tbb::parallel_for(0, mesh.vertexNum, 1, [&](int ii) {
            q[ii] += mesh.masses[ii] * c[ii];
            }
        );


        tbb::parallel_for(0, (int)BH.D1Index.size(), 1, [&](int ii) {
            Matrix3d H = BH.H3x3[ii];

            Vector3d tempC;

            //for (int i = 0; i < 3; i++) {
            tempC = c[BH.D1Index[ii]];
            //}
            Vector3d tempQ = H * tempC;


            countMutex[BH.D1Index[ii]].lock();
            q[BH.D1Index[ii]] += tempQ;
            countMutex[BH.D1Index[ii]].unlock();

            }
        );

        double tempSum = 0;

        tempSum = parallel_reduce(
            tbb::blocked_range<int>(0, mesh.vertexNum), 0.0,
            [&](const tbb::blocked_range<int>& rg, double temp_sum) {
                for (int i = rg.begin(); i != rg.end(); i++) {
                    q[i] = mesh.Constraints[i] * q[i];
                    temp_sum += (c[i].dot(q[i]));
                }
                return temp_sum;
            },
            [&](double left, double right) {
                return left + right;
            }
        );

        double alpha = deltaN / tempSum;
        //cout << "tempSum:------------------"<<tempSum << endl;
        //if(tempSum)
        deltaO = deltaN;
        deltaN = 0;
        vector<Vector3d> s(mesh.vertexNum, Vector3d(0, 0, 0));


        deltaN = parallel_reduce(
            tbb::blocked_range<int>(0, mesh.vertexNum), 0.0,
            [&](const tbb::blocked_range<int>& rg, double temp_deltaN) {
                for (int i = rg.begin(); i != rg.end(); i++) {
                    dX[i] = dX[i] + alpha * c[i];


                    r[i] = r[i] - alpha * q[i];


                    s[i] = P[i] * r[i];

                    temp_deltaN += (r[i].dot(s[i]));
                }
                return temp_deltaN;
            },
            [&](double left, double right) {
                return left + right;
            }
        );

        tbb::parallel_for(0, mesh.vertexNum, 1, [&](int i) {
            c[i] = s[i] + (deltaN / deltaO) * c[i];
            c[i] = mesh.Constraints[i] * c[i];
            }
        );

    }
    //printf("cg counts: %d\n", cgCounts);
    if (cgCounts == 0) {
        printf("indefinite exit\n");
        exit(0);
    }
    //if (localOptimal)
    //    return tempDeltaX;
    return dX;
}

void
calculateMovingDirection(const mesh3D& mesh, BHessian& BH, vector<Vector3d> gradient, vector<Vector3d>& direction) {
    direction = cholmod_solver(BH, gradient, mesh);
    //direction = Eigen_CG_solver(BH, gradient, mesh);
    //direction = PCG_Solver(mesh, BH, gradient);
}

void computeEnergyVal(const mesh3D& mesh, double& energyVal, const Ground& gd, double Kappa) {
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
                for (int j = 0; j != 4; ++j)
                    for (int k = 0; k != 4; ++k)
                        for (int d = 0; d != 3; ++d)
                            temp_bendingEnergyVal += verts[info.verts[j]][d] * info.Q(j, k) * verts[info.verts[k]][d];
            }
            return temp_bendingEnergyVal;
        },
        [&](double left, double right) {
            return left + right;
        });

    bendingEnergyVal *= (0.5 * mesh.bendingStiffness);
    energyVal += bendingEnergyVal;
#else
    energyVal += getObjBending_Energy(mesh);
#endif

    //energyVal += getObjEnergy_AniostroI5_3D(mesh.vertexes, mesh, lengthRate, contract_ratio);
    energyVal *= mesh.IPC_dt * mesh.IPC_dt;
    double deltaE = 0;

    deltaE = parallel_reduce(
        tbb::blocked_range<int>(0, mesh.vertexNum), 0.0,
        [&](const tbb::blocked_range<int>& rg, double temp_deltaE) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                temp_deltaE += ((mesh.vertexes[i] - mesh.xTilta[i]).squaredNorm() * mesh.masses[i] / 2.0);
                temp_deltaE += ((mesh.vertexes[i] - mesh.V_prev[i]).squaredNorm() * (mesh.drag_coeff) * mesh.masses[i] / 2.0);
            }
            return temp_deltaE;
        },
        [&](double left, double right) {
            return left + right;
        }
    );

    energyVal += deltaE;
    Eigen::VectorXd constraintVals, bVals;
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
    //TODO: parallelize
    tbb::parallel_for(startCI, (int)constraintVals.size(), 1, [&](int cI)
        //for (int cI = startCI; cI < constraintVals.size(); ++cI)
        {

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

    tbb::parallel_for(startCI, (int)constraintVals.size(), 1, [&](int cI)
        //for (int cI = startCI; cI < constraintVals.size(); ++cI)
        {

            const MMCVID& MMCVIDI = mesh.Self_EE_ActiveSet[cI - startCI];
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
    //double bE
    energyVal += Kappa * bVals.sum();
#ifdef USE_FRICTION
    // self friction
    double fricDHat = mesh.Fhat * mesh.IPC_dt * mesh.IPC_dt;
    double eps = std::sqrt(fricDHat);

    Eigen::VectorXd EI(mesh.Self_activeSet_lastH.size());

    tbb::parallel_for(0, (int)mesh.Self_activeSet_lastH.size(), 1, [&](int cI)
        {
            Eigen::RowVector3d relDX3D;

            const auto& MMCVIDI = mesh.Self_activeSet_lastH[cI];
            if (MMCVIDI[0] >= 0) {
                // edge-edge
                IPC::computeRelDX_EE(mesh.vertexes[MMCVIDI[0]] - mesh.V_prev[MMCVIDI[0]],
                    mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]],
                    mesh.vertexes[MMCVIDI[2]] - mesh.V_prev[MMCVIDI[2]],
                    mesh.vertexes[MMCVIDI[3]] - mesh.V_prev[MMCVIDI[3]],
                    mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], relDX3D);
            }
            else {
                // point-triangle and degenerate edge-edge
                assert(MMCVIDI[1] >= 0);
                if (MMCVIDI[2] < 0) {
                    // PP
                    IPC::computeRelDX_PP(mesh.vertexes[-MMCVIDI[0] - 1] - mesh.V_prev[-MMCVIDI[0] - 1],
                        mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]], relDX3D);
                }
                else if (MMCVIDI[3] < 0) {
                    // PE
                    IPC::computeRelDX_PE(mesh.vertexes[-MMCVIDI[0] - 1] - mesh.V_prev[-MMCVIDI[0] - 1],
                        mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]],
                        mesh.vertexes[MMCVIDI[2]] - mesh.V_prev[MMCVIDI[2]],
                        mesh.MMDistCoord[cI][0], relDX3D);
                }
                else {
                    // PT
                    IPC::computeRelDX_PT(mesh.vertexes[-MMCVIDI[0] - 1] - mesh.V_prev[-MMCVIDI[0] - 1],
                        mesh.vertexes[MMCVIDI[1]] - mesh.V_prev[MMCVIDI[1]],
                        mesh.vertexes[MMCVIDI[2]] - mesh.V_prev[MMCVIDI[2]],
                        mesh.vertexes[MMCVIDI[3]] - mesh.V_prev[MMCVIDI[3]],
                        mesh.MMDistCoord[cI][0], mesh.MMDistCoord[cI][1], relDX3D);
                }
            }

            EI[cI] = 0.0;
            double relDXSqNorm = (relDX3D * mesh.MMTanBasis[cI]).squaredNorm();
            if (relDXSqNorm > fricDHat) {
                EI[cI] += mesh.Self_lambda_lastH[cI] * std::sqrt(relDXSqNorm);
            }
            else {
                double f0;
                IPC::f0_SF(relDXSqNorm, eps, f0);
                EI[cI] += mesh.Self_lambda_lastH[cI] * f0;
            }
        }
    );

    energyVal += EI.sum() * mesh.friction;

    // ground friction
    double Ef = 0;
    int contactPairI = 0;
    for (const auto& vI : mesh.Environment_activeSet_lastH) {
        Eigen::Matrix<double, 3, 1> VDiff = (mesh.vertexes[vI] - mesh.V_prev[vI]).transpose();
        //        VDiff -= Base::velocitydt;
        Eigen::Matrix<double, 3, 1> VProj = VDiff - VDiff.dot(gd.normal) * gd.normal;
        double VProjMag2 = VProj.squaredNorm();
        if (VProjMag2 > fricDHat) {
            Ef += mesh.friction * mesh.Environment_lambda_lastH[contactPairI] * (std::sqrt(VProjMag2) - eps * 0.5);
        }
        else {
            Ef += mesh.friction * mesh.Environment_lambda_lastH[contactPairI] * VProjMag2 / eps * 0.5;
        }
        ++contactPairI;
    }
    //Ef *= 1.0;
    energyVal += Ef;
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
        std::unordered_set<int> sEdgeInds;
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
            sh.build(mesh, mesh.averageEdgeLength);
        }

        sh.calculateActivateSet(mesh);
    }
    gd.calculateActivateSet(mesh);
}

bool lineSearch(mesh3D& mesh,
    SpatialHash& sh,
    const Ground& gd,
    const vector<Vector3d>& searchDir,
    const vector<Vector3d>& gradient,
    double& stepSize,
    double armijoParam,
    double lowerBound,
    double Kappa) {
    //std::stringstream msg;



    bool stopped = false;
    double lastEnergyVal;

    computeEnergyVal(mesh, lastEnergyVal, gd, Kappa);

    //msg << "E_last = " << lastEnergyVal << " stepSize: " << stepSize << " -> ";

    // const double m = searchDir.dot(gradient);
    // const double c1m = 1.0e-4 * m;
    double c1m = 0.0;
    armijoParam = 0;
    if (armijoParam > 0.0) {
        c1m = parallel_reduce(
            tbb::blocked_range<int>(0, mesh.vertexNum), 0.0,
            [&](const tbb::blocked_range<int>& rg, double temp_deltaE) {
                for (int i = rg.begin(); i != rg.end(); i++) {
                    temp_deltaE += gradient[i].dot(searchDir[i]);
                }
                return temp_deltaE;
            },
            [&](double left, double right) {
                return left + right;
            }
        );
    }
    c1m *= armijoParam;
    vector<Vector3d> resultV0 = mesh.vertexes;

    stepForward(resultV0, searchDir, mesh, stepSize);


    int numOfIntersect = 0;
    if (mesh.use_barrier) {
        sh.build(mesh, mesh.averageEdgeLength);
        while (isIntersected(gd, sh, mesh, resultV0)) {
            numOfIntersect++;
            printf("intersect\n");
            stepSize /= 2.0;
            stepForward(resultV0, searchDir, mesh, stepSize);
            sh.build(mesh, mesh.averageEdgeLength);
        }
        //msg << stepSize << "(safeGuard_IP) -> ";





    }

    buildCollisionSets(mesh, sh, gd, false);
    double testingE;
    computeEnergyVal(mesh, testingE, gd, Kappa);


    //cout << "test test test:  " << abs(lastEnergyVal - mesh.restSNKE)<<"        "<< c1m<< endl;

    int numOfLineSearch = 0;
    double LFStepSize = stepSize;
    while ((testingE > lastEnergyVal + stepSize * c1m) && (stepSize > SimulationParameters::lineSearchMinFraction *
        LFStepSize)) {

        printf("ls iteration id:  %d,  testingE:  %f       lastEnergyVal:   %f\n", numOfLineSearch, testingE, lastEnergyVal);
        stepSize /= 2.0;
        //}

        ++numOfLineSearch;
        if (stepSize == 0.0) {
            stopped = true;
            break;
        }

        stepForward(resultV0, searchDir, mesh, stepSize);

        //msg << stepSize << " -> ";

        buildCollisionSets(mesh, sh, gd);
        computeEnergyVal(mesh, testingE, gd, Kappa);
    }

    if (stepSize < LFStepSize) {
        bool needRecomputeCS = false;
        if (mesh.use_barrier) {
            while (isIntersected(gd, sh, mesh, resultV0)) {
                stepSize /= 2.0;
                numOfIntersect++;
                printf("intersect\n");
                stepForward(resultV0, searchDir, mesh, stepSize);
                sh.build(mesh, mesh.averageEdgeLength);
                needRecomputeCS = true;
            }
        }
        if (needRecomputeCS) {
            buildCollisionSets(mesh, sh, gd, false);
        }
    }

    lastEnergyVal = testingE;

    return stopped;
}

void suggestKappa(double& kappa, const double& Hhat, const double& bboxDiagSize2, const double& meanMass) {
    double H_b;
    compute_H_b(1.0e-16 * bboxDiagSize2, Hhat, H_b);
    kappa = SimulationParameters::kappaMaxMassFactor * meanMass / (4.0e-16 * bboxDiagSize2 * H_b);
}

void upperBoundKappa(double& kappa, const double& Hhat, const double& bboxDiagSize2, const double& meanMass) {
    double H_b;
    //double bboxDiagSize2 = (maxCorner - minCorner).squaredNorm();
    compute_H_b(1.0e-16 * bboxDiagSize2, Hhat, H_b);
    double kappaMax = SimulationParameters::kappaUpperBoundFactor * SimulationParameters::kappaMaxMassFactor * meanMass / (4.0e-16 * bboxDiagSize2 * H_b);
    if (kappa > kappaMax) {
        kappa = kappaMax;
    }
    //printf("max Kappa: %f\n", kappaMax);
}

void initKappa(mesh3D& mesh, const Ground& grd, double& kappa) {
    std::vector<int> constraintStartInds;
    buildConstraintStartIndsWithMM(mesh.Environment_ActiveSet, mesh.Self_ActiveSet, constraintStartInds);

    if (constraintStartInds.back()) {
        vector<Vector3d> g_E(mesh.vertexNum, Vector3d(0, 0, 0)), g_c(mesh.vertexNum, Vector3d(0, 0, 0));

        computeEGradient(mesh, g_E);
        //compute_fiction_gradient(mesh, grd, g_E, mesh.Hhat * IPC_dt * IPC_dt, FEM::friction);
        VectorXd constraintVal;
        int startCI = constraintStartInds[0];
        Evaluate_GroundConstraintVals(grd, mesh, constraintVal, startCI);
        //animConfig.collisionObjects[coI]->evaluateConstraints(result, activeSet[coI], constraintVal);



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
#ifndef NEWB
        compute_g_dpt(mesh, mesh.Self_ActiveSet, constraintVal, g_c, startCI, 1);
#else
        compute_g_dpt_new(mesh, mesh.Self_ActiveSet, constraintVal, g_c, startCI, 1, mesh.Hhat);
#endif
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


void postLineSearch(mesh3D& mesh, const Ground& grd, double alpha, double& kappa) {
    if (kappa == 0.0) {
        initKappa(mesh, grd, kappa);
    }
    else {
        //TODO: avoid recomputation of constraint functions
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

        //double bboxDiagSize2 = (mesh.maxCorner - mesh.minCorner).squaredNorm();
        //double dTol = 1e-18 * mesh.bboxDiagSize2;
        //animConfig.collisionObjects[coI]->evaluateConstraints(result,
        //    activeSet[coI], constraintVal);
        //double bboxDiagSize2 = (mesh.maxCorner - mesh.minCorner).squaredNorm();
        for (int i = 0; i < mesh.Environment_ActiveSet.size(); ++i) {
            if (constraintVal[constraintValIndStart + i] < mesh.dTol) {
                mesh.closeConstraintID.emplace_back(mesh.Environment_ActiveSet[i]);
                mesh.closeConstraintVal.emplace_back(constraintVal[constraintValIndStart + i]);
            }
        }

        constraintValIndStart = constraintVal.size();
        Evaluate_SelfPTConstraintVals(mesh, constraintVal, constraintValIndStart);

        for (int i = 0; i < mesh.Self_ActiveSet.size(); ++i) {
            // std::cout << MMActiveSet.back()[i][0] << " " << MMActiveSet.back()[i][1] << " " << MMActiveSet.back()[i][2] << " " << MMActiveSet.back()[i][3] << ", d=" << constraintVal[constraintValIndStart + i] << std::endl;
            if (constraintVal[constraintValIndStart + i] < mesh.dTol) {
                mesh.closeMConstraintID.emplace_back(mesh.Self_ActiveSet[i]);
                mesh.closeMConstraintVal.emplace_back(constraintVal[constraintValIndStart + i]);
            }
        }
    }


}

void updateBoundaryMoveDir(mesh3D& mesh, vector<Vector3d>& moveDir, double ipc_dt, double alpha) {
    double angleX = PI / 2.5 * ipc_dt * alpha;
    Matrix3d rotationL, rotationR;
    rotationL << 1, 0, 0, 0, cos(angleX), sin(angleX), 0, -sin(angleX), cos(angleX);
    rotationR << 1, 0, 0, 0, cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);

    double mvl = -1 * ipc_dt * alpha;
    tbb::parallel_for(0, (int)(moveDir.size()), 1, [&](int i)
        //for (int i = 0; i < moveDir.size();i++)
        {
            if (mesh.boundaryTypes[i] == 1) {
                //moveDir[i] = rotationL * mesh.vertexes[i] - mesh.vertexes[i];
                moveDir[i] = Vector3d(mvl, 0, 0);

            }
            //                          if (mesh.boundaryTypes[i] == -1) {
            //                              moveDir[i] = rotationR * mesh.vertexes[i] - mesh.vertexes[i];
            //                          }
        }
    );
}


double calculate_distToOpt_PN(const vector<Vector3d>& moveDir) {
    double distToOpt_PN = 0;
    distToOpt_PN = parallel_reduce(
        tbb::blocked_range<int>(0, (int)moveDir.size()), 0.0,
        [&](const tbb::blocked_range<int>& rg, double temp_max) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                for (int jj = 0; jj < 3; jj++) {
                    if (temp_max < abs(moveDir[i][jj])) {
                        temp_max = abs(moveDir[i][jj]);
                    }
                }
            }
            return temp_max;
        },
        [&](double left, double right) {
            return left > right ? left : right;
        }
    );

    return distToOpt_PN;
}

static bool hasConverged(const mesh3D& mesh, const vector<Vector3d>& moveDir) {
    double distToOpt_PN = calculate_distToOpt_PN(moveDir);
    bool gradVanish = (distToOpt_PN < sqrt(SimulationParameters::newtonConvergenceTolFactor * mesh.bboxDiagSize2 * mesh.IPC_dt * mesh.IPC_dt));
    return gradVanish;
}

static bool computeNewtonDirection(mesh3D& mesh, BHessian& BH, vector<Vector3d>& gradient, vector<Vector3d>& moveDir) {
    calculateMovingDirection(mesh, BH, gradient, moveDir);
    return true;
}

static double computeFeasibleStepSize(mesh3D& mesh, SpatialHash& sh, Ground& gd, vector<Vector3d>& moveDir) {
    double alpha = 1.0;
    double slackness_a = SimulationParameters::environmentSlackness;
    double slackness_m = SimulationParameters::selfContactSlackness;

    //filterStepSize(mesh, moveDir, alpha);

    Environment_largestFeasibleStepSize(mesh, gd, moveDir, slackness_a, alpha);

    if (alpha <= 0) {
        alpha = 1;
    }
    std::vector<std::pair<int, int>> newCandidates;

    if (mesh.use_barrier) {

        Self_largestFeasibleStepSize(mesh, sh, moveDir, slackness_m, newCandidates, alpha);

        double partialCCD_alpha = alpha;

        Eigen::VectorXd pMag(mesh.surfVerts.size());
        for (int i = 0; i < pMag.size(); ++i) {
            int surfId = mesh.surfVerts[i];
            pMag[i] = moveDir[surfId].norm();
        }
        double alpha_CFL = std::sqrt(mesh.Hhat) / (pMag.maxCoeff() * 2.0);

        double fullCCD_alpha = alpha;
        sh.build(mesh, moveDir, fullCCD_alpha, mesh.averageEdgeLength);
        Self_largestFeasibleStepSize_CCD(mesh, sh, moveDir, slackness_m, fullCCD_alpha);

        alpha = min(alpha, alpha_CFL);

        if (partialCCD_alpha > 2 * alpha_CFL) {
            alpha = min(partialCCD_alpha, fullCCD_alpha * 1.0);
            alpha = max(alpha, alpha_CFL);
        }
    }

    return alpha;
}

static bool performLineSearch(mesh3D& mesh, SpatialHash& sh, Ground& gd, vector<Vector3d>& moveDir, vector<Vector3d>& gradient, double& alpha, double Kappa) {
    double alpha_feasible = alpha;

    bool isStop = lineSearch(mesh, sh, gd, moveDir, gradient, alpha, 0, 0, Kappa);
    return isStop;
}

static void finalizeNewtonStep(mesh3D& mesh, double alpha, double& totalTimeStep, vector<Vector3d>& moveDir) {
    totalTimeStep += alpha;
}

int solve_subIP(mesh3D& mesh, SpatialHash& sh, Ground& gd, double Kappa, float& time0, float& time1, float& time2, float& time3, float& time4, double& collisionNum) {
    int iterCap = SimulationParameters::newtonIterCap, k = 0;

    vector<Vector3d> moveDir(mesh.vertexNum, Vector3d(0, 0, 0));

    double totalTimeStep = 0;
    for (; k < iterCap; ++k) {

        HighResolutionTimerForWin timer0, timer1, timer2, timer3, timer4;

        vector<Vector3d> gradient(mesh.vertexNum, Vector3d(0, 0, 0));
        BHessian BH;
#if defined USE_QUADRATIC_BENDING
        size_t fem12x12Size = mesh.tetrahedraNum + mesh.quadBendingInfo.size();
#else
        size_t fem12x12Size = mesh.tetrahedraNum + mesh.tri_edges.size();
#endif
        BH.H3x3.reserve(mesh.Environment_ActiveSet.size() + mesh.Self_ActiveSet.size() * 4);
        BH.H6x6.reserve(mesh.Self_ActiveSet.size() * 2 + mesh.Self_EE_ActiveSet.size() * 2);
        BH.H9x9.reserve(mesh.Self_ActiveSet.size() + mesh.Self_EE_ActiveSet.size());
        BH.H12x12.reserve(fem12x12Size + mesh.Self_ActiveSet.size() + mesh.Self_EE_ActiveSet.size() * 2);
        BH.D1Index.reserve(mesh.Environment_ActiveSet.size());
        BH.D2Index.reserve(mesh.Self_EE_ActiveSet.size());
        BH.D3Index.reserve(mesh.Self_ActiveSet.size());
        BH.D4Index.reserve(mesh.Self_ActiveSet.size() + mesh.Self_EE_ActiveSet.size());
        timer0.set_start();
        collisionNum += computeGradientAndHessian(mesh, gradient, BH, gd);

        timer0.set_end();
        timer1.set_start();

        if (k && hasConverged(mesh, moveDir)/* && totalTimeStep > 1 - 1e-3*/) {
            break;
        }

        computeNewtonDirection(mesh, BH, gradient, moveDir);

        timer1.set_end();
        timer2.set_start();

        double alpha = computeFeasibleStepSize(mesh, sh, gd, moveDir);

        timer2.set_end();
        timer3.set_start();

        bool isStop = performLineSearch(mesh, sh, gd, moveDir, gradient, alpha, Kappa);
        timer3.set_end();
        timer4.set_start();
        postLineSearch(mesh, gd, alpha, Kappa);
        timer4.set_end();

        float time00, time11, time22, time33, time44;
        time00 = timer0.get_millisecond();
        time11 = timer1.get_millisecond();
        time22 = timer2.get_millisecond();
        time33 = timer3.get_millisecond();
        time44 = timer4.get_millisecond();

        time0 += time00;
        time1 += time11;
        time2 += time22;
        time3 += time33;
        time4 += time44;
        finalizeNewtonStep(mesh, alpha, totalTimeStep, moveDir);
        //printf("time0 = %f,  time1 = %f,  time2 = %f,  time3 = %f,  time4 = %f\n", time00, time11, time22, time33, time44);
    }
    printf("newton iteration:  %d    and    Kappa:  %f\n", k, mesh.Kappa);
    total_iter += k;
    return k;
}

void export_obj(const mesh3D& mesh, int index) {
    std::ofstream cloth_stream("output/surface_obj_" + std::to_string(index) + ".obj");
    cloth_stream << "# Generated by hc"
        << "\n";
    cloth_stream << std::fixed << std::setprecision(6) << "cloth\n";
    auto cloth_vertices = mesh.vertexes;
    auto cloth_triangles = mesh.surface;

    for (auto vec : cloth_vertices) {
        cloth_stream << "v " << vec(0) << " " << vec(1) << " "
            << vec(2) << "\n";
    }
    cloth_stream << "s 1\n";
    for (auto tri : cloth_triangles) {
        cloth_stream << "f " << tri(0) + 1 << " " << tri(1) + 1
            << " " << tri(2) + 1 << "\n";
    }
    cloth_stream.close();
}

bool isRotate = false;
bool loadTempTimeInfo = true;

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


int IPC_Solver(int& stepId, mesh3D& mesh, SpatialHash& sh, Ground& gd) {
    //mesh3D& mesh = meshTetes->mesh3Ds[0];

    if (loadTempTimeInfo) {
        std::string fileVertex = "tempData/timeCost.txt";
        ifstream ifs(fileVertex);
        if (ifs) {
            ifs >> ttime0 >> ttime1 >> ttime2 >> ttime3 >> ttime4 >> time_total >> total_iter >> step_index
                >> mesh.Kappa;
            loadTempTimeInfo = false;
        }
        ifs.close();
    }

    stepId = step_index;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //printf("hello 0 p\n");
    upperBoundKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
    if (mesh.Kappa < 1e-16) {
        suggestKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
    }
    initKappa(mesh, gd, mesh.Kappa);
    //printf("hello 1 p\n");
#ifdef USE_FRICTION
    buildFriction(mesh, gd);
#endif

    if (true && mesh.update_hard_constraint_functor != nullptr) {
        vector<Vector3d> moveDir(mesh.vertexNum, Vector3d(0, 0, 0));
        double new_alpha = 1;
        //for (int i = 0; i < mesh.boundary_vertexes_indices.size(); i++)
        int boundary_vertex_num = mesh.boundary_vertexes_indices.size();
        tbb::parallel_for(0, boundary_vertex_num, 1, [&](int i)
            {
                moveDir[mesh.boundary_vertexes_indices[i]] = mesh.update_hard_constraint_functor(
                    mesh.vertexes[mesh.boundary_vertexes_indices[i]], new_alpha, mesh.IPC_dt);
            }
        );
        sh.build(mesh, moveDir, new_alpha, mesh.averageEdgeLength);
        Self_largestFeasibleStepSize_CCD(mesh, sh, moveDir, 0.8, new_alpha);
        tbb::parallel_for(0, boundary_vertex_num, 1, [&](int i)
            {
                moveDir[mesh.boundary_vertexes_indices[i]] = mesh.update_hard_constraint_functor(
                    mesh.vertexes[mesh.boundary_vertexes_indices[i]], new_alpha, mesh.IPC_dt);
            }
        );
        vector<Vector3d> resultV0 = mesh.vertexes;
        stepForward(resultV0, moveDir, mesh, 1, true);

        sh.build(mesh, mesh.averageEdgeLength);
        int numOfIntersect = 0;

        while (isIntersected(gd, sh, mesh, mesh.vertexes)) {
            new_alpha /= 2.0;
            tbb::parallel_for(0, boundary_vertex_num, 1, [&](int i)
                {
                    moveDir[mesh.boundary_vertexes_indices[i]] = mesh.update_hard_constraint_functor(
                        mesh.vertexes[mesh.boundary_vertexes_indices[i]], new_alpha, mesh.IPC_dt);
                }
            );
            stepForward(resultV0, moveDir, mesh, 1, true);
            sh.build(mesh, mesh.averageEdgeLength);
        }
        cout << "new_alpha:";
        cout << new_alpha << endl;
        buildCollisionSets(mesh, sh, gd, false);
    }
    //printf("rotate\n");
    int k = 0;
    float time0 = 0;
    float time1 = 0;
    float time2 = 0;
    float time3 = 0;
    float time4 = 0;
    double collisonNum = 0;
    while (true) {
        mesh.closeConstraintID.resize(0);
        mesh.closeMConstraintID.resize(0);
        mesh.closeConstraintVal.resize(0);
        mesh.closeMConstraintVal.resize(0);
        k = solve_subIP(mesh, sh, gd, mesh.Kappa, time0, time1, time2, time3, time4, collisonNum);

        VectorXd constraintVals;
        int offset = 0;
        Evaluate_GroundConstraintVals(gd, mesh, constraintVals, offset);
        offset = constraintVals.size();
        Evaluate_SelfPTConstraintVals(mesh, constraintVals, offset);

        if (constraintVals.size()) {
            double minm = constraintVals.minCoeff();
            double maxm = constraintVals.maxCoeff();
            //std::cout << minm << "    " << maxm << endl;
            if (constraintVals.minCoeff() < mesh.dTol) {
                break;
            }
            else if (constraintVals.maxCoeff() < mesh.Hhat) {
                break;
            }
        }
        else {
            break;
        }
    }


    updateVelocity(mesh.vertexes, mesh.V_prev, mesh.velocities, mesh.IPC_dt, mesh.vertexNum, mesh.is_quasi_static);


    mesh.V_prev = mesh.vertexes;
    computeXTilta(mesh);

    ttime0 += time0;
    ttime1 += time1;
    ttime2 += time2;
    ttime3 += time3;
    ttime4 += time4;

    std::cout << "                                                finished a step" << endl;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time for a step = "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    time_total += std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    totalCollision += collisonNum;

    std::cout << "totalCollision  = "
        << totalCollision << std::endl;


    std::cout << "averageCollision  = "
        << totalCollision / total_iter << std::endl;


    std::cout << "total  = "
        << time_total << "[ms]" << std::endl;
    std::cout << "total newton step:" << total_iter << std::endl;

    std::cout << "frame number  = "
        << step_index++ << std::endl;

    ofstream outTime("timeCost.txt");

    outTime << "time0: " << ttime0 / 1000.0 << endl;
    outTime << "time1: " << ttime1 / 1000.0 << endl;
    outTime << "time2: " << ttime2 / 1000.0 << endl;
    outTime << "time3: " << ttime3 / 1000.0 << endl;
    outTime << "time4: " << ttime4 / 1000.0 << endl;

    outTime << "timeAll: " << time_total / 1000.0 << endl;
    outTime << "total iter: " << total_iter << endl;
    outTime << "frames: " << step_index << endl;
    outTime << "totalCollisionNum: " << totalCollision << endl;

    outTime << "averageCollision: " << totalCollision / total_iter << endl;
    outTime.close();

    if (step_index % 10 == 0)
    {
        ofstream outTime2("tempData/timeCost.txt");
        outTime2 << ttime0 << endl;
        outTime2 << ttime1 << endl;
        outTime2 << ttime2 << endl;
        outTime2 << ttime3 << endl;
        outTime2 << ttime4 << endl;
        outTime2 << time_total << endl;
        outTime2 << total_iter << endl;
        outTime2 << step_index << endl;
        outTime2 << mesh.Kappa << endl;
        outTime2.close();
        mesh.output_tetTempData();
    }

    //if (step_index == 250) {
    //    exit(0);
    //}
    return k;
}