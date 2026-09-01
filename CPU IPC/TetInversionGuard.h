#pragma once

#include "SimulationMesh.h"

#include <vector>

// Limits stepSize so no tetrahedron shrinks below volumeRatio times its
// current signed volume under x <- x - alpha * searchDirection.
void limitStepToPreventTetInversion(
    const mesh3D& mesh,
    const std::vector<Eigen::Vector3d>& searchDirection,
    double volumeRatio,
    double& stepSize);
