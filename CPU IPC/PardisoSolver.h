#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <array>
#include <cstddef>
#include <vector>

class PardisoSolver final {
public:
    PardisoSolver();
    ~PardisoSolver();

    PardisoSolver(const PardisoSolver&) = delete;
    PardisoSolver& operator=(const PardisoSolver&) = delete;

    void solve(
        const Eigen::SparseMatrix<double>& lowerTriangle,
        const Eigen::VectorXd& rightHandSide,
        Eigen::VectorXd& result,
        int threadCount);

    std::size_t symbolicAnalysisCount() const { return symbolicAnalysisCount_; }
    std::size_t numericFactorizationCount() const { return numericFactorizationCount_; }
    int factorNonZeros() const { return factorNonZeros_; }
    int activeThreadCount() const { return activeThreadCount_; }
    double analysisMilliseconds() const { return analysisMilliseconds_; }
    double factorizationMilliseconds() const { return factorizationMilliseconds_; }
    double solveMilliseconds() const { return solveMilliseconds_; }

private:
    bool hasSamePattern(const Eigen::SparseMatrix<double>& matrix) const;
    void initialize();
    void release();
    void runPhase(
        int phase,
        const Eigen::SparseMatrix<double>& matrix,
        double* rightHandSide,
        double* result);
    void configureThreads(int threadCount);

    std::array<void*, 64> internalData_{};
    std::array<int, 64> parameters_{};
    std::vector<int> permutation_;
    std::vector<int> outerIndices_;
    std::vector<int> innerIndices_;
    Eigen::VectorXd rightHandSideBuffer_;
    int dimension_ = 0;
    int factorNonZeros_ = 0;
    int activeThreadCount_ = 0;
    int defaultThreadCount_ = 0;
    bool initialized_ = false;
    bool permutationReady_ = false;
    bool forceReordering_ = false;
    bool updateFactorReference_ = false;
    int permutationDimension_ = 0;
    int referenceFactorNonZeros_ = 0;
    double analysisMilliseconds_ = 0.0;
    double factorizationMilliseconds_ = 0.0;
    double solveMilliseconds_ = 0.0;
    std::size_t symbolicAnalysisCount_ = 0;
    std::size_t numericFactorizationCount_ = 0;
};
