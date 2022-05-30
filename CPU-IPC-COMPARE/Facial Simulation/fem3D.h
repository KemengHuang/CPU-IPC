#pragma once
#ifndef FEM3D_H
#define FEM3D_H
//#include "mesh.h"
#include "fem_math.h"
#include "collisionUtil.h"

//host functions
void initMesh3D(mesh3D& mesh, int type, double scale);

double calculateVolum(const vector<Vector3d>& vertexes, const Vector4i& index);
Matrix3d calculateDms3D_double(const vector<Vector3d>& vertexes, const Vector4i& index, const int& i);
double getObjEnergy_StableNHK2_3D(const vector<Vector3d>& vertexes, const mesh3D& mesh, const double& lengthRate, const double& volumeRate);
double getObjRestEnergy_StableNHK2_3D(const vector<Vector3d>& vertexes, const mesh3D& mesh, const double& lengthRate, const double& volumeRate);
MatrixXd computePFPX3D_double(const Matrix3d& InverseDm);
Matrix3d computePEPF_StableNHK3D_double(const Matrix3d& F, const double& lengthRate, const double& volumRate);
Matrix3d computePEPF_StableNHK3D_2_double(const Matrix3d& F, const double& lengthRate, const double& volumRate);
Matrix3d computePEPF_Aniostropic3D_double(const Matrix3d& F, Vector3d direction, const double& scale, const double& contract_length);
Matrix3d computePEPF_AniostropicRehabi3D_double(const Matrix3d& F, Vector3d direction, const double& u_anios);
MatrixXd project_StabbleNHK_H_3D(const Matrix3d& F, const double& lengthRate, const double& volumRate);
MatrixXd project_StabbleNHK_2_H_3D(const Matrix3d& F, const double& lengthRate, const double& volumRate);
MatrixXd project_ANIOSI5_H_3D(const Matrix3d& F, Vector3d direction, const double& scale, const double& contract_length);
MatrixXd project_ANIOSI5_Rehabi_H_3D(const Matrix3d& F, Vector3d direction, const double& u_anios);



#endif // !FEM3D_CUH



