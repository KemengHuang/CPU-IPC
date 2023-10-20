#pragma once
#ifndef FEM3D_H
#define FEM3D_H
//#include "mesh.h"
#include "fem_math.h"
#include "collisionUtil.h"

//host functions
void initMesh3D(mesh3D &mesh, int type, double scale);

double calculateVolum(const vector<Vector3d> &vertexes, const Vector4i &index);

void __Inverse(const Eigen::Matrix3d &input, Eigen::Matrix3d &result);

void __Inverse2x2(const Eigen::Matrix2d &input, Eigen::Matrix2d &result);

Matrix3d calculateDms3D_double(const vector<Vector3d> &vertexes, const Vector4i &index, const int &i);

Matrix<double, 3, 2> calculateDs32D_double(const vector<Vector3d> &vertexes, const Vector3i &index);

double getObjEnergy_StableNHK2_3D(const vector<Vector3d> &vertexes, const mesh3D &mesh, const double &lengthRate,
                                  const double &volumeRate);

double getObjRestEnergy_StableNHK2_3D(const vector<Vector3d> &vertexes, const mesh3D &mesh, const double &lengthRate,
                                      const double &volumeRate);

MatrixXd computePFPX3D_double(const Matrix3d &InverseDm);

MatrixXd computePFPX32D_double(const Matrix2d &InverseDm);

Matrix3d computePEPF_StableNHK3D_double(const Matrix3d &F, const double &lengthRate, const double &volumRate);

Matrix<double, 3, 2>
computePEPF_baraffwitkin_double(const Matrix<double, 3, 2> &F,
                                const Vector2d &anisotropic_a,
                                const Vector2d &anisotropic_b,
                                double stretchS, double shearS);


Matrix<double, 6, 6>
project_baraffwitkint_H_3D(const Matrix<double, 3, 2> &F,
                                const Vector2d &anisotropic_a,
                                const Vector2d &anisotropic_b,
                                double stretchS, double shearS);
double getObjEnergy_baraffwitkin_3D(const mesh3D &mesh,
                                  const Vector2d &anisotropic_a,
                                  const Vector2d &anisotropic_b,
                                  double stretchS, double shearS);

Matrix3d computePEPF_StableNHK3D_2_double(const Matrix3d &F, const double &lengthRate, const double &volumRate);

Matrix3d computePEPF_Aniostropic3D_double(const Matrix3d &F, Vector3d direction, const double &scale,
                                          const double &contract_length);

Matrix3d computePEPF_AniostropicRehabi3D_double(const Matrix3d &F, Vector3d direction, const double &u_anios);

MatrixXd project_StabbleNHK_H_3D(const Matrix3d &F, const double &lengthRate, const double &volumRate);

MatrixXd project_StabbleNHK_2_H_3D(const Matrix3d &F, const double &lengthRate, const double &volumRate);

MatrixXd
project_ANIOSI5_H_3D(const Matrix3d &F, Vector3d direction, const double &scale, const double &contract_length);

MatrixXd project_ANIOSI5_Rehabi_H_3D(const Matrix3d &F, Vector3d direction, const double &u_anios);


#endif // !FEM3D_CUH



