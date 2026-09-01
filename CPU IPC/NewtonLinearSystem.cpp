#include "NewtonLinearSystem.h"

#include "ContactMechanics.h"
#include "SimulationMesh.h"

#include <Eigen/IterativeLinearSolvers>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

NewtonLinearSystem::NewtonLinearSystem(int vertexCount)
    : matrix_(3 * vertexCount, 3 * vertexCount),
      rightHandSide_(3 * vertexCount),
      solution_(3 * vertexCount)
{
}

void NewtonLinearSystem::solve(
    const mesh3D& mesh,
    const BHessian& hessian,
    const std::vector<Eigen::Vector3d>& gradient,
    std::vector<Eigen::Vector3d>& direction,
    const LinearSolverOptions& options,
    std::size_t& matrixNonZeros)
{
    const std::size_t vertexCount = static_cast<std::size_t>(mesh.vertexNum);
    if (gradient.size() != vertexCount
        || direction.size() != vertexCount
        || mesh.boundaryTypes.size() != vertexCount
        || mesh.masses.size() != vertexCount
        || matrix_.rows() != 3 * mesh.vertexNum) {
        throw std::invalid_argument(
            "Newton linear-system buffers do not match the mesh vertex count");
    }

    assemble(mesh, hessian, gradient, matrixNonZeros);

    switch (options.backend) {
    case LinearSolverBackend::Cholmod:
        solveWithCholmod();
        break;
    case LinearSolverBackend::SuiteSparseLDL:
        solveWithSuiteSparseLDL();
        break;
    case LinearSolverBackend::EigenConjugateGradient:
        solveWithEigenConjugateGradient(options);
        break;
    default:
        throw std::invalid_argument("unsupported Newton linear solver backend");
    }

    scatter(direction);
}

std::size_t NewtonLinearSystem::symbolicAnalysisCount() const
{
    return cholmodSolver_.symbolicAnalysisCount()
        + suiteSparseLDLSolver_.symbolicAnalysisCount();
}

std::size_t NewtonLinearSystem::numericFactorizationCount() const
{
    return cholmodSolver_.numericFactorizationCount()
        + suiteSparseLDLSolver_.numericFactorizationCount();
}

void NewtonLinearSystem::assemble(
    const mesh3D& mesh,
    const BHessian& hessian,
    const std::vector<Eigen::Vector3d>& gradient,
    std::size_t& matrixNonZeros)
{
    hessian.toTriplets(mesh.boundaryTypes, triplets_);
    const int massOffset = static_cast<int>(triplets_.size());
    triplets_.resize(massOffset + 3 * mesh.vertexNum);

    tbb::parallel_for(0, 3 * mesh.vertexNum, 1, [&](int degreeOfFreedom) {
        triplets_[massOffset + degreeOfFreedom] = Eigen::Triplet<double>(
            degreeOfFreedom,
            degreeOfFreedom,
            mesh.masses[degreeOfFreedom / 3] * (1.0 + mesh.drag_coeff));
    });

    matrix_.setZero();
    matrix_.setFromTriplets(triplets_.begin(), triplets_.end());
    matrixNonZeros = (std::max)(
        matrixNonZeros, static_cast<std::size_t>(matrix_.nonZeros()));

    tbb::parallel_for(0, static_cast<int>(gradient.size()), 1, [&](int vertex) {
        rightHandSide_.template segment<3>(3 * vertex) =
            mesh.boundaryTypes[vertex] == 0
                ? gradient[vertex]
                : Eigen::Vector3d::Zero();
    });
}

void NewtonLinearSystem::solveWithCholmod()
{
    cholmodSolver_.set_pattern(matrix_);
    cholmodSolver_.solve(rightHandSide_, solution_);
}

void NewtonLinearSystem::solveWithSuiteSparseLDL()
{
    suiteSparseLDLSolver_.solve(matrix_, rightHandSide_, solution_);
}

void NewtonLinearSystem::solveWithEigenConjugateGradient(
    const LinearSolverOptions& options)
{
    if (!std::isfinite(options.relativeTolerance)
        || options.relativeTolerance <= 0.0
        || options.maximumIterations <= 0) {
        throw std::invalid_argument(
            "Eigen conjugate-gradient tolerances must be positive");
    }

    using Preconditioner = Eigen::IncompleteCholesky<double, Eigen::Lower>;
    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<double>,
        Eigen::Lower,
        Preconditioner> solver;
    solver.setTolerance(options.relativeTolerance);
    solver.setMaxIterations(options.maximumIterations);
    solver.compute(matrix_);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error(
            "Eigen conjugate-gradient preconditioner setup failed");
    }

    solution_ = solver.solve(rightHandSide_);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error(
            "Eigen conjugate-gradient solve failed after "
            + std::to_string(solver.iterations())
            + " iterations; estimated error="
            + std::to_string(solver.error()));
    }
}

void NewtonLinearSystem::scatter(
    std::vector<Eigen::Vector3d>& direction) const
{
    tbb::parallel_for(0, static_cast<int>(direction.size()), 1, [&](int vertex) {
        direction[vertex] = solution_.template segment<3>(3 * vertex);
    });
}
