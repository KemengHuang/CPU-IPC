#pragma once
#ifndef FEM3D_H
#define FEM3D_H
//#include "mesh/Mesh.h"
#include "math/MathUtils.h"
#include "collision/SpatialHash.h"


double edgeTheta(
    const Eigen::Vector3d& q0,
    const Eigen::Vector3d& q1,
    const Eigen::Vector3d& q2,
    const Eigen::Vector3d& q3,
    Eigen::Matrix<double, 1, 12>* derivative, // edgeVertex, then edgeOppositeVertex
    Eigen::Matrix<double, 12, 12>* hessian);

//host functions
void initMesh3D(mesh3D &mesh, int type, double scale);

double calculateVolum(const std::vector<Eigen::Vector3d> &vertexes, const Eigen::Vector4i &index);

Eigen::Matrix3d calculateDms3D_double(const std::vector<Eigen::Vector3d> &vertexes, const Eigen::Vector4i &index, const int &i);

Eigen::Matrix<double, 3, 2> calculateDs32D_double(const std::vector<Eigen::Vector3d> &vertexes, const Eigen::Vector3i &index);

double getObjEnergy_StableNHK2_3D(const std::vector<Eigen::Vector3d> &vertexes, const mesh3D &mesh, const double &lengthRate,
                                  const double &volumeRate);

double getObjRestEnergy_StableNHK2_3D(const std::vector<Eigen::Vector3d> &vertexes, const mesh3D &mesh, const double &lengthRate,
                                      const double &volumeRate);

Eigen::MatrixXd computePFPX3D_double(const Eigen::Matrix3d &InverseDm);

Eigen::MatrixXd computePFPX32D_double(const Eigen::Matrix2d &InverseDm);

Eigen::Matrix3d computePEPF_StableNHK3D_double(const Eigen::Matrix3d &F, const double &lengthRate, const double &volumRate);

Eigen::Matrix<double, 3, 2>
computePEPF_baraffwitkin_double(const Eigen::Matrix<double, 3, 2>& F,
    const Eigen::Vector2d& anisotropic_a,
    const Eigen::Vector2d& anisotropic_b,
    double stretchS, double shearS, double strainRate);


Eigen::Matrix<double, 6, 6>
project_baraffwitkint_H_3D(const Eigen::Matrix<double, 3, 2>& F,
    const Eigen::Vector2d& anisotropic_a,
    const Eigen::Vector2d& anisotropic_b,
    double stretchS, double shearS, double strainRate);
double getObjEnergy_baraffwitkin_3D(const mesh3D& mesh,
    const Eigen::Vector2d& anisotropic_a,
    const Eigen::Vector2d& anisotropic_b,
    double stretchS, double shearS, double strainRate);

double getObjBending_Energy(const mesh3D& mesh);

Eigen::Matrix3d computePEPF_StableNHK3D_2_double(const Eigen::Matrix3d &F, const double &lengthRate, const double &volumRate);

Eigen::Matrix3d computePEPF_Aniostropic3D_double(const Eigen::Matrix3d &F, Eigen::Vector3d direction, const double &scale,
                                          const double &contract_length);

Eigen::Matrix3d computePEPF_AniostropicRehabi3D_double(const Eigen::Matrix3d &F, Eigen::Vector3d direction, const double &u_anios);

Eigen::MatrixXd project_StabbleNHK_H_3D(const Eigen::Matrix3d &F, const double &lengthRate, const double &volumRate);

Eigen::MatrixXd project_StabbleNHK_2_H_3D(const Eigen::Matrix3d &F, const double &lengthRate, const double &volumRate);

Eigen::MatrixXd
project_ANIOSI5_H_3D(const Eigen::Matrix3d &F, Eigen::Vector3d direction, const double &scale, const double &contract_length);

Eigen::MatrixXd project_ANIOSI5_Rehabi_H_3D(const Eigen::Matrix3d &F, Eigen::Vector3d direction, const double &u_anios);


#endif // !FEM3D_CUH



