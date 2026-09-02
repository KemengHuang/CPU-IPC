#pragma once

#include "LinearSolverOptions.h"
#include "CholmodSolver.h"
#ifdef CIPC_HAS_PARDISO
#include "PardisoSolver.h"
#endif

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cstddef>
#include <vector>

class BHessian;
class mesh3D;

class NewtonLinearSystem final {
public:
    explicit NewtonLinearSystem(int vertexCount);

    void solve(
        const mesh3D& mesh,
        const BHessian& hessian,
        const std::vector<Eigen::Vector3d>& gradient,
        std::vector<Eigen::Vector3d>& direction,
        const LinearSolverOptions& options,
        std::size_t& matrixNonZeros);

    std::size_t symbolicAnalysisCount() const;
    std::size_t numericFactorizationCount() const;
    double pardisoAnalysisMilliseconds() const;
    double pardisoFactorizationMilliseconds() const;
    double pardisoSolveMilliseconds() const;
    int cholmodThreadCount() const;
    int pardisoThreadCount() const;
    int pardisoFactorNonZeros() const;
    int vertexCount() const { return static_cast<int>(matrix_.rows() / 3); }

private:
    void assemble(
        const mesh3D& mesh,
        const BHessian& hessian,
        const std::vector<Eigen::Vector3d>& gradient,
        std::size_t& matrixNonZeros);
    void solveWithCholmod(const LinearSolverOptions& options);
    void solveWithPardiso(const LinearSolverOptions& options);
    void solveWithEigenConjugateGradient(const LinearSolverOptions& options);
    void scatter(std::vector<Eigen::Vector3d>& direction) const;

    std::vector<Eigen::Triplet<double>> triplets_;
    Eigen::SparseMatrix<double> matrix_;
    Eigen::VectorXd rightHandSide_;
    Eigen::VectorXd solution_;
    CholmodSolver cholmodSolver_;
#ifdef CIPC_HAS_PARDISO
    PardisoSolver pardisoSolver_;
#endif
};
