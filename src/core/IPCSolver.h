#pragma once
#ifndef _IPC_FUNC_H_
#define _IPC_FUNC_H_

#include "fem/Elasticity.h"

// Global solver state.  Originally these were loose top-level variables in
// IPCSolver.cpp; grouping them makes the state explicit and simplifies future
// work (e.g. checkpointing, multi-solver instances).
struct IPCSolverState {
    int step_index = 0;
    double time_total = 0.0;
    double ttime0 = 0.0;
    double ttime1 = 0.0;
    double ttime2 = 0.0;
    double ttime3 = 0.0;
    double ttime4 = 0.0;
    double totalCollision = 0.0;
    int total_iter = 0;
};

extern IPCSolverState g_solverState;

int solve_subIP(mesh3D& mesh, SpatialHash& sh, Ground& gd, double Kappa, float &time0, float& time1, float& time2, float& time3, float& time4, double& collisionNum);
int IPC_Solver(int& stepId, mesh3D& mesh, SpatialHash& sh, Ground& gd);
void computeXTilta(mesh3D& mesh);
#endif // !_IPC_FUNC_H_
