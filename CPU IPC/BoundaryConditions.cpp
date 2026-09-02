#include "BoundaryConditions.h"

#include "ContactMechanics.h"
#include "SimulationMesh.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace {

void validateVertexIndices(
    const std::vector<int>& indices,
    std::size_t vertexCount,
    const char* label)
{
    std::unordered_set<int> uniqueIndices;
    uniqueIndices.reserve(indices.size());
    for (int vertex : indices) {
        if (vertex < 0 || static_cast<std::size_t>(vertex) >= vertexCount) {
            throw std::invalid_argument(
                std::string(label) + " contains an out-of-range vertex index");
        }
        if (!uniqueIndices.insert(vertex).second) {
            throw std::invalid_argument(
                std::string(label) + " contains a duplicate vertex index");
        }
    }
}

void validateSoftState(const mesh3D& mesh)
{
    const SoftBoundaryCondition& soft = mesh.boundaryConditions.soft;
    if (soft.vertexIndices.empty()) {
        if (!soft.targetPositions.empty() || soft.updateTarget) {
            throw std::invalid_argument(
                "soft boundary targets and updates require soft boundary vertices");
        }
        return;
    }
    if (!std::isfinite(soft.weight) || soft.weight <= 0.0) {
        throw std::invalid_argument("soft boundary weight must be finite and positive");
    }
    if (soft.targetPositions.size() != soft.vertexIndices.size()) {
        throw std::invalid_argument(
            "soft boundary targets must match the soft vertex count");
    }
    for (const Eigen::Vector3d& target : soft.targetPositions) {
        if (!target.allFinite()) {
            throw std::invalid_argument("soft boundary target must be finite");
        }
    }
}

} // namespace

namespace BoundaryConditionOps {

void initialize(mesh3D& mesh)
{
    const std::size_t vertexCount = mesh.vertexes.size();
    if (mesh.vertexNum < 0
        || static_cast<std::size_t>(mesh.vertexNum) != vertexCount
        || mesh.v_rest.size() != vertexCount
        || mesh.boundaryTypes.size() != vertexCount) {
        throw std::invalid_argument(
            "boundary conditions require complete vertex, rest, and type arrays");
    }

    for (int type : mesh.boundaryTypes) {
        if (type < boundaryTypeCode(VertexBoundaryType::Free)) {
            throw std::invalid_argument("boundary type codes must be non-negative");
        }
    }

    DirichletBoundaryCondition& dirichlet = mesh.boundaryConditions.dirichlet;
    validateVertexIndices(dirichlet.vertexIndices, vertexCount, "Dirichlet boundary");
    if (static_cast<bool>(dirichlet.updateDirection)
        != !dirichlet.vertexIndices.empty()) {
        throw std::invalid_argument(
            "animated Dirichlet vertices and their update functor must be configured together");
    }
    for (int vertex : dirichlet.vertexIndices) {
        if (isFreeBoundary(mesh.boundaryTypes[vertex])) {
            throw std::invalid_argument(
                "animated Dirichlet vertices must be marked constrained");
        }
    }

    SoftBoundaryCondition& soft = mesh.boundaryConditions.soft;
    validateVertexIndices(soft.vertexIndices, vertexCount, "soft boundary");
    for (int vertex : soft.vertexIndices) {
        if (!isFreeBoundary(mesh.boundaryTypes[vertex])) {
            throw std::invalid_argument(
                "soft boundary vertices must remain free Newton unknowns");
        }
    }
    if (!soft.vertexIndices.empty() && soft.targetPositions.empty()) {
        soft.targetPositions.resize(soft.vertexIndices.size());
        for (std::size_t index = 0; index < soft.vertexIndices.size(); ++index) {
            soft.targetPositions[index] = mesh.vertexes[soft.vertexIndices[index]];
        }
    }
    validateSoftState(mesh);
}

void updateSoftTargets(mesh3D& mesh, int stepIndex)
{
    SoftBoundaryCondition& soft = mesh.boundaryConditions.soft;
    if (soft.vertexIndices.empty() || !soft.updateTarget) {
        return;
    }
    if (stepIndex < 0 || !std::isfinite(mesh.IPC_dt) || mesh.IPC_dt <= 0.0) {
        throw std::invalid_argument("soft boundary update received invalid time data");
    }
    validateSoftState(mesh);

    for (std::size_t index = 0; index < soft.vertexIndices.size(); ++index) {
        const int vertex = soft.vertexIndices[index];
        const Eigen::Vector3d target = soft.updateTarget(
            mesh.vertexes[vertex],
            mesh.v_rest[vertex],
            stepIndex,
            mesh.IPC_dt);
        if (!target.allFinite()) {
            throw std::runtime_error("soft boundary update produced a non-finite target");
        }
        soft.targetPositions[index] = target;
    }
}

void evaluateDirichletDirections(
    const mesh3D& mesh,
    int stepIndex,
    double stepFraction,
    std::vector<Eigen::Vector3d>& directions)
{
    const DirichletBoundaryCondition& dirichlet =
        mesh.boundaryConditions.dirichlet;
    if (!dirichlet.updateDirection || dirichlet.vertexIndices.empty()) {
        throw std::logic_error("animated Dirichlet boundary is not configured");
    }
    if (stepIndex < 0
        || !std::isfinite(stepFraction)
        || stepFraction < 0.0 || stepFraction > 1.0
        || !std::isfinite(mesh.IPC_dt) || mesh.IPC_dt <= 0.0) {
        throw std::invalid_argument("Dirichlet update received invalid time data");
    }

    directions.assign(mesh.vertexes.size(), Eigen::Vector3d::Zero());
    tbb::parallel_for(
        0, static_cast<int>(dirichlet.vertexIndices.size()), 1, [&](int index) {
            const int vertex = dirichlet.vertexIndices[index];
            directions[vertex] = dirichlet.updateDirection(
                mesh.vertexes[vertex],
                mesh.v_rest[vertex],
                stepIndex,
                stepFraction,
                mesh.IPC_dt);
        });

    for (int vertex : dirichlet.vertexIndices) {
        if (!directions[vertex].allFinite()) {
            throw std::runtime_error(
                "Dirichlet boundary update produced a non-finite direction");
        }
    }
}

void addSoftDerivatives(
    const mesh3D& mesh,
    std::vector<Eigen::Vector3d>& gradient,
    BHessian& hessian)
{
    const SoftBoundaryCondition& soft = mesh.boundaryConditions.soft;
    if (soft.vertexIndices.empty()) {
        return;
    }
    validateSoftState(mesh);
    if (gradient.size() != mesh.vertexes.size()) {
        throw std::invalid_argument(
            "soft boundary gradient does not match the mesh vertex count");
    }
    if (hessian.H3x3.size() != hessian.D1Index.size()) {
        throw std::logic_error("soft boundary received inconsistent 3x3 Hessian storage");
    }

    const std::size_t hessianOffset = hessian.H3x3.size();
    hessian.H3x3.resize(hessianOffset + soft.vertexIndices.size());
    hessian.D1Index.resize(hessianOffset + soft.vertexIndices.size());
    const Eigen::Matrix3d localHessian =
        soft.weight * Eigen::Matrix3d::Identity();

    tbb::parallel_for(
        0, static_cast<int>(soft.vertexIndices.size()), 1, [&](int index) {
            const int vertex = soft.vertexIndices[index];
            gradient[vertex] += soft.weight
                * (mesh.vertexes[vertex] - soft.targetPositions[index]);
            hessian.H3x3[hessianOffset + index] = localHessian;
            hessian.D1Index[hessianOffset + index] = vertex;
        });
}

double softEnergy(const mesh3D& mesh)
{
    const SoftBoundaryCondition& soft = mesh.boundaryConditions.soft;
    if (soft.vertexIndices.empty()) {
        return 0.0;
    }
    validateSoftState(mesh);

    return tbb::parallel_reduce(
        tbb::blocked_range<int>(0, static_cast<int>(soft.vertexIndices.size())),
        0.0,
        [&](const tbb::blocked_range<int>& range, double energy) {
            for (int index = range.begin(); index != range.end(); ++index) {
                const int vertex = soft.vertexIndices[index];
                energy += 0.5 * soft.weight
                    * (mesh.vertexes[vertex] - soft.targetPositions[index]).squaredNorm();
            }
            return energy;
        },
        [](double left, double right) { return left + right; });
}

} // namespace BoundaryConditionOps
