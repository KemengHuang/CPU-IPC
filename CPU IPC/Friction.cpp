#include "Friction.h"
#include "FrictionKinematics.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <tbb/parallel_for.h>

namespace {

template <int Size, typename LiftToMesh>
Eigen::Matrix<double, Size, Size> analyticFrictionHessian(
    const Eigen::Vector2d& relativeDisplacement,
    double epsilon,
    double multiplier,
    LiftToMesh liftToMesh)
{
    using VectorN = Eigen::Matrix<double, Size, 1>;
    using MatrixN = Eigen::Matrix<double, Size, Size>;

    MatrixN hessian = MatrixN::Zero();
    const double weight = (std::max)(0.0, multiplier);
    if (weight == 0.0) {
        return hessian;
    }

    auto addEigenDirection = [&](const Eigen::Vector2d& tangentDirection, double eigenvalue) {
        const double nonnegativeEigenvalue = (std::max)(0.0, eigenvalue);
        if (nonnegativeEigenvalue == 0.0) {
            return;
        }
        VectorN lifted;
        liftToMesh(tangentDirection, lifted);
        hessian.noalias() += weight * nonnegativeEigenvalue
            * lifted * lifted.transpose();
    };

    const double squaredNorm = relativeDisplacement.squaredNorm();
    const double norm = std::sqrt(squaredNorm);
    if (norm > 0.0) {
        const Eigen::Vector2d parallel = relativeDisplacement / norm;
        const Eigen::Vector2d perpendicular(-parallel.y(), parallel.x());
        if (norm > epsilon) {
            // Hessian of ||u||: zero along u and 1/||u|| orthogonal to u.
            addEigenDirection(perpendicular, 1.0 / norm);
        }
        else {
            double perpendicularEigenvalue = 0.0;
            double parallelEigenvalue = 0.0;
            IPC::f1_SF_div_relDXNorm(
                squaredNorm, epsilon, perpendicularEigenvalue);
            IPC::f2_SF(squaredNorm, epsilon, parallelEigenvalue);
            addEigenDirection(parallel, parallelEigenvalue);
            addEigenDirection(perpendicular, perpendicularEigenvalue);
        }
    }
    else {
        double isotropicEigenvalue = 0.0;
        IPC::f1_SF_div_relDXNorm(0.0, epsilon, isotropicEigenvalue);
        addEigenDirection(Eigen::Vector2d::UnitX(), isotropicEigenvalue);
        addEigenDirection(Eigen::Vector2d::UnitY(), isotropicEigenvalue);
    }
    return hessian;
}

Eigen::Matrix<double, 3, 2> groundTangentBasis(const Eigen::Vector3d& normal)
{
    const Eigen::Vector3d helper = std::abs(normal.x()) < 0.9
        ? Eigen::Vector3d::UnitX()
        : Eigen::Vector3d::UnitY();
    Eigen::Matrix<double, 3, 2> basis;
    basis.col(0) = (helper - helper.dot(normal) * normal).normalized();
    basis.col(1) = normal.cross(basis.col(0)).normalized();
    return basis;
}

} // namespace

namespace Friction {

void initialize(mesh3D& mesh, const Ground& grd) {
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


void addGradient(mesh3D& mesh, const Ground& grd, std::vector<Eigen::Vector3d>& grad_inc, double eps2,
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

void addHessian(
    mesh3D& mesh,
    const Ground& grd,
    BHessian& BH,
    double eps2,
    double coef)
{
    const double eps = std::sqrt(eps2);

    // Ground friction uses the C0 clamp implemented by its energy/gradient.
    const Eigen::Matrix<double, 3, 2> groundBasis =
        groundTangentBasis(grd.normal);
    for (int contact = 0;
         contact < static_cast<int>(mesh.Environment_activeSet_lastH.size());
         ++contact) {
        const int vertex = mesh.Environment_activeSet_lastH[contact];
        const double multiplier = (std::max)(
            0.0, coef * mesh.Environment_lambda_lastH[contact]);
        const Eigen::Vector2d relativeDisplacement =
            groundBasis.transpose() * (mesh.vertexes[vertex] - mesh.V_prev[vertex]);
        const double norm = relativeDisplacement.norm();

        Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
        if (multiplier > 0.0) {
            if (norm > eps) {
                const Eigen::Vector2d parallel = relativeDisplacement / norm;
                const Eigen::Vector2d perpendicular(-parallel.y(), parallel.x());
                const Eigen::Vector3d lifted = groundBasis * perpendicular;
                hessian.noalias() = (multiplier / norm)
                    * lifted * lifted.transpose();
            }
            else {
                hessian.noalias() = (multiplier / eps)
                    * (groundBasis.col(0) * groundBasis.col(0).transpose()
                        + groundBasis.col(1) * groundBasis.col(1).transpose());
            }
        }
        BH.H3x3.emplace_back(hessian);
        BH.D1Index.emplace_back(vertex);
    }

    // For frozen lambda, closest-point coordinates, and tangent basis, the
    // self-friction Hessian has an analytic two-dimensional eigensystem.
    for (int contact = 0;
         contact < static_cast<int>(mesh.Self_activeSet_lastH.size());
         ++contact) {
        const EncodedContact& ids = mesh.Self_activeSet_lastH[contact];
        const double multiplier = coef * mesh.Self_lambda_lastH[contact];
        Eigen::RowVector3d relativeDisplacement3D;

        if (ids[0] >= 0) {
            IPC::computeRelDX_EE(
                mesh.vertexes[ids[0]] - mesh.V_prev[ids[0]],
                mesh.vertexes[ids[1]] - mesh.V_prev[ids[1]],
                mesh.vertexes[ids[2]] - mesh.V_prev[ids[2]],
                mesh.vertexes[ids[3]] - mesh.V_prev[ids[3]],
                mesh.MMDistCoord[contact][0],
                mesh.MMDistCoord[contact][1],
                relativeDisplacement3D);
            const Eigen::Vector2d relativeDisplacement =
                (relativeDisplacement3D * mesh.MMTanBasis[contact]).transpose();
            const auto hessian = analyticFrictionHessian<12>(
                relativeDisplacement,
                eps,
                multiplier,
                [&](const Eigen::Vector2d& tangentDirection,
                    Eigen::Matrix<double, 12, 1>& lifted) {
                    IPC::liftRelDXTanToMesh_EE(
                        tangentDirection,
                        mesh.MMTanBasis[contact],
                        mesh.MMDistCoord[contact][0],
                        mesh.MMDistCoord[contact][1],
                        lifted);
                });
            BH.H12x12.emplace_back(hessian);
            BH.D4Index.emplace_back(ids[0], ids[1], ids[2], ids[3]);
            continue;
        }

        const int point = -ids[0] - 1;
        if (ids[2] < 0) {
            IPC::computeRelDX_PP(
                mesh.vertexes[point] - mesh.V_prev[point],
                mesh.vertexes[ids[1]] - mesh.V_prev[ids[1]],
                relativeDisplacement3D);
            const Eigen::Vector2d relativeDisplacement =
                (relativeDisplacement3D * mesh.MMTanBasis[contact]).transpose();
            const auto hessian = analyticFrictionHessian<6>(
                relativeDisplacement,
                eps,
                multiplier,
                [&](const Eigen::Vector2d& tangentDirection,
                    Eigen::Matrix<double, 6, 1>& lifted) {
                    IPC::liftRelDXTanToMesh_PP(
                        tangentDirection, mesh.MMTanBasis[contact], lifted);
                });
            BH.H6x6.emplace_back(hessian);
            BH.D2Index.emplace_back(point, ids[1]);
        }
        else if (ids[3] < 0) {
            IPC::computeRelDX_PE(
                mesh.vertexes[point] - mesh.V_prev[point],
                mesh.vertexes[ids[1]] - mesh.V_prev[ids[1]],
                mesh.vertexes[ids[2]] - mesh.V_prev[ids[2]],
                mesh.MMDistCoord[contact][0],
                relativeDisplacement3D);
            const Eigen::Vector2d relativeDisplacement =
                (relativeDisplacement3D * mesh.MMTanBasis[contact]).transpose();
            const auto hessian = analyticFrictionHessian<9>(
                relativeDisplacement,
                eps,
                multiplier,
                [&](const Eigen::Vector2d& tangentDirection,
                    Eigen::Matrix<double, 9, 1>& lifted) {
                    IPC::liftRelDXTanToMesh_PE(
                        tangentDirection,
                        mesh.MMTanBasis[contact],
                        mesh.MMDistCoord[contact][0],
                        lifted);
                });
            BH.H9x9.emplace_back(hessian);
            BH.D3Index.emplace_back(point, ids[1], ids[2]);
        }
        else {
            IPC::computeRelDX_PT(
                mesh.vertexes[point] - mesh.V_prev[point],
                mesh.vertexes[ids[1]] - mesh.V_prev[ids[1]],
                mesh.vertexes[ids[2]] - mesh.V_prev[ids[2]],
                mesh.vertexes[ids[3]] - mesh.V_prev[ids[3]],
                mesh.MMDistCoord[contact][0],
                mesh.MMDistCoord[contact][1],
                relativeDisplacement3D);
            const Eigen::Vector2d relativeDisplacement =
                (relativeDisplacement3D * mesh.MMTanBasis[contact]).transpose();
            const auto hessian = analyticFrictionHessian<12>(
                relativeDisplacement,
                eps,
                multiplier,
                [&](const Eigen::Vector2d& tangentDirection,
                    Eigen::Matrix<double, 12, 1>& lifted) {
                    IPC::liftRelDXTanToMesh_PT(
                        tangentDirection,
                        mesh.MMTanBasis[contact],
                        mesh.MMDistCoord[contact][0],
                        mesh.MMDistCoord[contact][1],
                        lifted);
                });
            BH.H12x12.emplace_back(hessian);
            BH.D4Index.emplace_back(point, ids[1], ids[2], ids[3]);
        }
    }
}

double energy(
    const mesh3D& mesh,
    const Ground& ground,
    EnergyWorkspace& workspace)
{
    const double epsilonSquared = mesh.Fhat * mesh.IPC_dt * mesh.IPC_dt;
    const double epsilon = std::sqrt(epsilonSquared);

    Eigen::VectorXd& contactEnergies = workspace.contactEnergies;
    contactEnergies.resize(mesh.Self_activeSet_lastH.size());
    tbb::parallel_for(
        0,
        static_cast<int>(mesh.Self_activeSet_lastH.size()),
        1,
        [&](int contact) {
            const EncodedContact& ids = mesh.Self_activeSet_lastH[contact];
            Eigen::RowVector3d relativeDisplacement3D;
            if (ids[0] >= 0) {
                IPC::computeRelDX_EE(
                    mesh.vertexes[ids[0]] - mesh.V_prev[ids[0]],
                    mesh.vertexes[ids[1]] - mesh.V_prev[ids[1]],
                    mesh.vertexes[ids[2]] - mesh.V_prev[ids[2]],
                    mesh.vertexes[ids[3]] - mesh.V_prev[ids[3]],
                    mesh.MMDistCoord[contact][0],
                    mesh.MMDistCoord[contact][1],
                    relativeDisplacement3D);
            }
            else {
                const int point = -ids[0] - 1;
                if (ids[2] < 0) {
                    IPC::computeRelDX_PP(
                        mesh.vertexes[point] - mesh.V_prev[point],
                        mesh.vertexes[ids[1]] - mesh.V_prev[ids[1]],
                        relativeDisplacement3D);
                }
                else if (ids[3] < 0) {
                    IPC::computeRelDX_PE(
                        mesh.vertexes[point] - mesh.V_prev[point],
                        mesh.vertexes[ids[1]] - mesh.V_prev[ids[1]],
                        mesh.vertexes[ids[2]] - mesh.V_prev[ids[2]],
                        mesh.MMDistCoord[contact][0],
                        relativeDisplacement3D);
                }
                else {
                    IPC::computeRelDX_PT(
                        mesh.vertexes[point] - mesh.V_prev[point],
                        mesh.vertexes[ids[1]] - mesh.V_prev[ids[1]],
                        mesh.vertexes[ids[2]] - mesh.V_prev[ids[2]],
                        mesh.vertexes[ids[3]] - mesh.V_prev[ids[3]],
                        mesh.MMDistCoord[contact][0],
                        mesh.MMDistCoord[contact][1],
                        relativeDisplacement3D);
                }
            }

            const double squaredNorm =
                (relativeDisplacement3D * mesh.MMTanBasis[contact]).squaredNorm();
            if (squaredNorm > epsilonSquared) {
                contactEnergies[contact] =
                    mesh.Self_lambda_lastH[contact] * std::sqrt(squaredNorm);
            }
            else {
                double clampedNorm = 0.0;
                IPC::f0_SF(squaredNorm, epsilon, clampedNorm);
                contactEnergies[contact] =
                    mesh.Self_lambda_lastH[contact] * clampedNorm;
            }
        });

    double result = mesh.friction * contactEnergies.sum();
    for (int contact = 0;
         contact < static_cast<int>(mesh.Environment_activeSet_lastH.size());
         ++contact) {
        const int vertex = mesh.Environment_activeSet_lastH[contact];
        const Eigen::Vector3d difference =
            mesh.vertexes[vertex] - mesh.V_prev[vertex];
        const Eigen::Vector3d projected =
            difference - difference.dot(ground.normal) * ground.normal;
        const double squaredNorm = projected.squaredNorm();
        const double clampedNorm = squaredNorm > epsilonSquared
            ? std::sqrt(squaredNorm) - 0.5 * epsilon
            : 0.5 * squaredNorm / epsilon;
        result += mesh.friction
            * mesh.Environment_lambda_lastH[contact] * clampedNorm;
    }
    return result;
}

} // namespace Friction
