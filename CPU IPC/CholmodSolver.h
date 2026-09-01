#pragma once

#ifndef CIPC_CHOLMOD_SOLVER_H
#define CIPC_CHOLMOD_SOLVER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cholmod.h>
#include <cstddef>

class CholmodSolver final {
public:
    CholmodSolver();
    ~CholmodSolver();

    CholmodSolver(const CholmodSolver&) = delete;
    CholmodSolver& operator=(const CholmodSolver&) = delete;

    // Update numeric values and invalidate symbolic analysis only when the
    // compressed sparse structure actually changes.
    void set_pattern(const Eigen::SparseMatrix<double>& matrix);

    void solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& result);

    void preFactorize(const Eigen::SparseMatrix<double>& matrix);
    void solve_with_preFactorize(const Eigen::VectorXd& rhs, Eigen::VectorXd& result);

    std::size_t symbolicAnalysisCount() const { return symbolicAnalysisCount_; }
    std::size_t numericFactorizationCount() const { return numericFactorizationCount_; }

private:
    bool hasSamePattern(const Eigen::SparseMatrix<double>& matrix) const;
    void updateMatrixView();
    void factorize();
    void solveFactorized(const Eigen::VectorXd& rhs, Eigen::VectorXd& result);

    int numRows_ = 0;
    int numCols_ = 0;
    Eigen::VectorXi outerIndices_;
    Eigen::VectorXi innerIndices_;
    Eigen::VectorXd values_;

    cholmod_common common_;
    cholmod_sparse matrix_;
    cholmod_factor* factor_ = nullptr;
    cholmod_dense* solutionBuffer_ = nullptr;
    cholmod_dense* solveWorkspaceY_ = nullptr;
    cholmod_dense* solveWorkspaceE_ = nullptr;
    bool matrixReady_ = false;
    std::size_t symbolicAnalysisCount_ = 0;
    std::size_t numericFactorizationCount_ = 0;
};

#endif // CIPC_CHOLMOD_SOLVER_H
