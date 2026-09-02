#pragma once
#ifndef CIPC_BOUNDARY_CONDITIONS_H
#define CIPC_BOUNDARY_CONDITIONS_H

#include <Eigen/Core>

#include <cstddef>
#include <functional>
#include <vector>

class BHessian;
class mesh3D;

// boundaryTypes remains an integer array for compatibility with the mesh and
// collision kernels. These named values and predicates are the only places
// that should interpret its codes.
enum class VertexBoundaryType : int {
    Free = 0,
    Dirichlet = 1,
    ExternalCollider = 2
};

inline int boundaryTypeCode(VertexBoundaryType type)
{
    return static_cast<int>(type);
}

inline bool isFreeBoundary(int type)
{
    return type == boundaryTypeCode(VertexBoundaryType::Free);
}

inline bool isExternalColliderBoundary(int type)
{
    return type >= boundaryTypeCode(VertexBoundaryType::ExternalCollider);
}

struct DirichletBoundaryCondition {
    using UpdateFunctor = std::function<Eigen::Vector3d(
        const Eigen::Vector3d& currentPosition,
        const Eigen::Vector3d& restPosition,
        int stepIndex,
        double stepFraction,
        double timeStep)>;

    // Only animated Dirichlet vertices belong here. Static constrained
    // vertices need only a non-free entry in mesh3D::boundaryTypes. The
    // callback returns the solver direction p in x_new = x - p.
    std::vector<int> vertexIndices;
    UpdateFunctor updateDirection;
};

struct SoftBoundaryCondition {
    using UpdateFunctor = std::function<Eigen::Vector3d(
        const Eigen::Vector3d& currentPosition,
        const Eigen::Vector3d& restPosition,
        int stepIndex,
        double timeStep)>;

    // A target-position penalty: 0.5 * weight * ||x - target||^2. These
    // vertices remain free unknowns and therefore must have boundary type 0.
    std::vector<int> vertexIndices;
    std::vector<Eigen::Vector3d> targetPositions;
    double weight = 0.0;
    UpdateFunctor updateTarget;
};

struct BoundaryConditionSet {
    DirichletBoundaryCondition dirichlet;
    SoftBoundaryCondition soft;
};

namespace BoundaryConditionOps {

// Validate cross-array invariants and initialize omitted static soft targets
// from the current positions. Call once after mesh3D::v_rest is finalized.
void initialize(mesh3D& mesh);

// Refresh time-dependent soft targets once at the beginning of a time step.
void updateSoftTargets(mesh3D& mesh, int stepIndex);

// Evaluate the prescribed Dirichlet search direction at a fractional update.
// The output is always resized and cleared before constrained entries are set.
void evaluateDirichletDirections(
    const mesh3D& mesh,
    int stepIndex,
    double stepFraction,
    std::vector<Eigen::Vector3d>& directions);

// Add the mutually consistent gradient/Hessian of the soft penalty.
void addSoftDerivatives(
    const mesh3D& mesh,
    std::vector<Eigen::Vector3d>& gradient,
    BHessian& hessian);

double softEnergy(const mesh3D& mesh);

} // namespace BoundaryConditionOps

#endif // CIPC_BOUNDARY_CONDITIONS_H
