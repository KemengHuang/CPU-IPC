#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cstddef>
#include <vector>

class SuiteSparseLDLSolver final {
public:
    void solve(
        const Eigen::SparseMatrix<double>& lowerTriangle,
        const Eigen::VectorXd& rightHandSide,
        Eigen::VectorXd& result);

    std::size_t symbolicAnalysisCount() const { return symbolicAnalysisCount_; }
    std::size_t numericFactorizationCount() const { return numericFactorizationCount_; }

private:
    bool hasSamePattern(const Eigen::SparseMatrix<double>& lowerTriangle) const;
    void analyzePattern(const Eigen::SparseMatrix<double>& lowerTriangle);
    void updateValues(const Eigen::SparseMatrix<double>& lowerTriangle);
    void factorize();
    void solveFactorized(const Eigen::VectorXd& rightHandSide, Eigen::VectorXd& result);

    Eigen::SparseMatrix<double> matrix_;
    std::vector<int> outerIndices_;
    std::vector<int> innerIndices_;
    std::vector<int> columnPointers_;
    std::vector<int> eliminationTree_;
    std::vector<int> columnNonZeros_;
    std::vector<int> flags_;
    std::vector<int> permutation_;
    std::vector<int> inversePermutation_;
    std::vector<Eigen::Index> valueIndices_;
    std::vector<int> rowIndices_;
    std::vector<int> patternWorkspace_;
    std::vector<double> factorValues_;
    std::vector<double> diagonal_;
    std::vector<double> numericWorkspace_;
    bool patternReady_ = false;
    std::size_t symbolicAnalysisCount_ = 0;
    std::size_t numericFactorizationCount_ = 0;
};
