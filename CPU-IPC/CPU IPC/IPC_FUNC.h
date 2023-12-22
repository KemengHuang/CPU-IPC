#pragma once
#ifndef _IPC_FUNC_H_
#define _IPC_FUNC_H_
//#include <vector>
//#include "collisionUtil.h"
//#include "Eigen/Eigen"
#include "fem3D.h"

//void Evaluate_GroundConstraintVals(const Ground& gd, const mesh3D& mesh, Eigen::VectorXd& vals, const int& offset);
//void Evaluate_SelfConstraintVals(const mesh3D& mesh, Eigen::VectorXd& vals, const int& offset);
//double SelfConstraintVal(const mesh3D& mesh, const MMCVID& active);

int solve_subIP(mesh3D& mesh, SpatialHash& sh, Ground& gd, double Kappa, float &time0, float& time1, float& time2, float& time3, float& time4, double& collisionNum);
int IPC_Solver(int& stepId, mesh3D& mesh, SpatialHash& sh, Ground& gd);
void computeXTilta(mesh3D& mesh);
#endif // !_IPC_FUNC_H_

