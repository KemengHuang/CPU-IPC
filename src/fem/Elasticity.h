#pragma once
#ifndef FEM3D_H
#define FEM3D_H
#include "math/MathUtils.h"
#include "collision/SpatialHash.h"

double edgeTheta(
    const Eigen::Vector3d& q0,
    const Eigen::Vector3d& q1,
    const Eigen::Vector3d& q2,
    const Eigen::Vector3d& q3,
    Eigen::Matrix<double, 1, 12>* derivative, // edgeVertex, then edgeOppositeVertex
    Eigen::Matrix<double, 12, 12>* hessian);

// ------------------------------------------------------------------
// Mesh initialization
// ------------------------------------------------------------------
void initMesh3D(mesh3D& mesh, int type, double scale);

// ------------------------------------------------------------------
// Deformation gradient helpers
// ------------------------------------------------------------------
double calculateVolume(const std::vector<Eigen::Vector3d>& vertexes, const Eigen::Vector4i& index);

Eigen::Matrix3d calculateDms3D_double(const std::vector<Eigen::Vector3d>& vertexes, const Eigen::Vector4i& index, const int& i);

Eigen::Matrix<double, 3, 2> calculateDs32D_double(const std::vector<Eigen::Vector3d>& vertexes, const Eigen::Vector3i& index);

Eigen::Matrix<double, 9, 12> computePFPX3D_double(const Eigen::Matrix3d& InverseDm);

Eigen::Matrix<double, 6, 9> computePFPX32D_double(const Eigen::Matrix2d& InverseDm);

// ------------------------------------------------------------------
// Tetrahedral elasticity (Stable Neo-Hookean)
// ------------------------------------------------------------------
double getObjEnergy_StableNHK2_3D(const std::vector<Eigen::Vector3d>& vertexes, const mesh3D& mesh, const double& lengthRate,
                                  const double& volumeRate);

double getObjRestEnergy_StableNHK2_3D(const std::vector<Eigen::Vector3d>& vertexes, const mesh3D& mesh, const double& lengthRate,
                                      const double& volumeRate);

Eigen::Matrix3d computePEPF_StableNHK3D_double(const Eigen::Matrix3d& F, const double& lengthRate, const double& volumRate);

Eigen::Matrix3d computePEPF_StableNHK3D_2_double(const Eigen::Matrix3d& F, const double& lengthRate, const double& volumRate);

Eigen::MatrixXd project_StableNHK_H_3D(const Eigen::Matrix3d& F, const double& lengthRate, const double& volumRate);

Eigen::Matrix<double, 9, 9> project_StableNHK_2_H_3D(const Eigen::Matrix3d& F, const double& lengthRate, const double& volumRate);

// ------------------------------------------------------------------
// Triangle shell elasticity (Baraff-Witkin)
// ------------------------------------------------------------------
Eigen::Matrix<double, 3, 2>
computePEPF_baraffWitkin_double(const Eigen::Matrix<double, 3, 2>& F,
    const Eigen::Vector2d& anisotropic_a,
    const Eigen::Vector2d& anisotropic_b,
    double stretchS, double shearS, double strainRate);

Eigen::Matrix<double, 6, 6>
project_baraffWitkin_H_3D(const Eigen::Matrix<double, 3, 2>& F,
    const Eigen::Vector2d& anisotropic_a,
    const Eigen::Vector2d& anisotropic_b,
    double stretchS, double shearS, double strainRate);

double getObjEnergy_baraffWitkin_3D(const mesh3D& mesh,
    const Eigen::Vector2d& anisotropic_a,
    const Eigen::Vector2d& anisotropic_b,
    double stretchS, double shearS, double strainRate);

// ------------------------------------------------------------------
// Anisotropic muscle models
// ------------------------------------------------------------------
Eigen::Matrix3d computePEPF_Anisotropic3D_double(const Eigen::Matrix3d& F, const Eigen::Vector3d& direction, const double& scale,
                                                 const double& contract_length);

Eigen::Matrix3d computePEPF_AnisotropicRehabi3D_double(const Eigen::Matrix3d& F, const Eigen::Vector3d& direction, const double& u_anios);

Eigen::MatrixXd
project_ANIOSI5_H_3D(const Eigen::Matrix3d& F, const Eigen::Vector3d& direction, const double& scale, const double& contract_length);

Eigen::MatrixXd project_ANIOSI5_Rehabi_H_3D(const Eigen::Matrix3d& F, const Eigen::Vector3d& direction, const double& u_anios);

// ------------------------------------------------------------------
// Bending
// ------------------------------------------------------------------
double getObjBending_Energy(const mesh3D& mesh);

#endif // !FEM3D_CUH
