#pragma once
#ifndef _IPC_FUNC_H_
#define _IPC_FUNC_H_

#include "fem/Elasticity.h"

int solve_subIP(mesh3D& mesh, SpatialHash& sh, Ground& gd, double Kappa, float &time0, float& time1, float& time2, float& time3, float& time4, double& collisionNum);
int IPC_Solver(int& stepId, mesh3D& mesh, SpatialHash& sh, Ground& gd);
void computeXTilta(mesh3D& mesh);
#endif // !_IPC_FUNC_H_

