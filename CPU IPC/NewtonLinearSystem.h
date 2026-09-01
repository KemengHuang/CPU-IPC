#pragma once

#include "LinearSolverOptions.h"
#include "CholmodSolver.h"
#include "SuiteSparseLDLSolver.h"

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

private:
    void assemble(
        const mesh3D& mesh,
        const BHessian& hessian,
        const std::vector<Eigen::Vector3d>& gradient,
        std::size_t& matrixNonZeros);
    void solveWithCholmod();
    void solveWithSuiteSparseLDL();
    void solveWithEigenConjugateGradient(const LinearSolverOptions& options);
    void scatter(std::vector<Eigen::Vector3d>& direction) const;

    std::vector<Eigen::Triplet<double>> triplets_;
    Eigen::SparseMatrix<double> matrix_;
    Eigen::VectorXd rightHandSide_;
    Eigen::VectorXd solution_;
    CholmodSolver cholmodSolver_;
    SuiteSparseLDLSolver suiteSparseLDLSolver_;
};
