#pragma once
#ifndef CIPC_FEASIBLE_STEP_H
#define CIPC_FEASIBLE_STEP_H

#include "SimulationMesh.h"
#include "CollisionBroadPhase.h"

void limitStepByGround(
    const mesh3D& mesh,
    const Ground& ground,
    const std::vector<Eigen::Vector3d>& searchDirection,
    double slackness,
    double& stepSize);

#endif // CIPC_FEASIBLE_STEP_H
