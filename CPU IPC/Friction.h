#pragma once

#include "ContactMechanics.h"

namespace Friction {

struct EnergyWorkspace {
    Eigen::VectorXd contactEnergies;
};

void initialize(mesh3D& mesh, const Ground& ground);
void addGradient(
    mesh3D& mesh,
    const Ground& ground,
    std::vector<Eigen::Vector3d>& gradient,
    double epsilonSquared,
    double frictionCoefficient);
void addHessian(
    mesh3D& mesh,
    const Ground& ground,
    BHessian& hessian,
    double epsilonSquared,
    double frictionCoefficient);
double energy(
    const mesh3D& mesh,
    const Ground& ground,
    EnergyWorkspace& workspace);

} // namespace Friction
