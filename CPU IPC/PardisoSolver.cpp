#include "PardisoSolver.h"

#include <mkl.h>
#include <tbb/global_control.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace {

constexpr int realSymmetricPositiveDefinite = 2;
constexpr int analysisPhase = 11;
constexpr int factorizationPhase = 22;
constexpr int solvePhase = 33;
constexpr int releasePhase = -1;

template <typename Scalar>
bool arraysEqual(const Scalar* lhs, const Scalar* rhs, Eigen::Index count)
{
    return count == 0
        || std::memcmp(lhs, rhs, static_cast<std::size_t>(count) * sizeof(Scalar)) == 0;
}

} // namespace

PardisoSolver::PardisoSolver()
{
    static_assert(std::is_same<MKL_INT, int>::value,
        "PardisoSolver requires the oneMKL LP64 interface");
    static_assert(
        std::is_same<Eigen::SparseMatrix<double>::StorageIndex, int>::value,
        "PardisoSolver requires 32-bit Eigen sparse indices");
    initialize();
    defaultThreadCount_ = mkl_get_max_threads();
}

PardisoSolver::~PardisoSolver()
{
    release();
}

void PardisoSolver::initialize()
{
    internalData_.fill(nullptr);
    parameters_.fill(0);
    const int matrixType = realSymmetricPositiveDefinite;
    pardisoinit(internalData_.data(), &matrixType, parameters_.data());

    parameters_[0] = 1;   // Use the explicitly configured parameters below.
    parameters_[1] = 2;   // METIS nested-dissection ordering.
    parameters_[3] = 0;   // Direct factorization, no iterative-direct mode.
    parameters_[4] = 0;   // Do not provide a user permutation.
    parameters_[5] = 0;   // Write the solution to x and preserve b.
    parameters_[7] = 2;   // At most two iterative-refinement steps.
    parameters_[9] = 13;  // Pivot perturbation threshold used by oneMKL defaults.
    parameters_[10] = 0;  // Disable nonsymmetric scaling for this SPD matrix.
    parameters_[12] = 0;  // Disable weighted matching for this SPD matrix.
    parameters_[17] = -1; // Report factor nonzeros.
    parameters_[18] = -1; // Report factorization operations.
    parameters_[26] = 0;  // Input is validated by the project assembly path.
    parameters_[27] = 0;  // Double precision.
    parameters_[34] = 1;  // Zero-based C indexing.
    parameters_[36] = 0;  // CSR storage.
    parameters_[59] = 0;  // In-core factorization.
}

void PardisoSolver::release()
{
    if (!initialized_) {
        return;
    }

    const int maxFactorizations = 1;
    const int matrixNumber = 1;
    const int matrixType = realSymmetricPositiveDefinite;
    const int rightHandSideCount = 1;
    const int messageLevel = 0;
    const int phase = releasePhase;
    int error = 0;
    pardiso(
        internalData_.data(),
        &maxFactorizations,
        &matrixNumber,
        &matrixType,
        &phase,
        &dimension_,
        nullptr,
        nullptr,
        nullptr,
        permutation_.data(),
        &rightHandSideCount,
        parameters_.data(),
        &messageLevel,
        nullptr,
        nullptr,
        &error);
    initialized_ = false;
    internalData_.fill(nullptr);
}

bool PardisoSolver::hasSamePattern(
    const Eigen::SparseMatrix<double>& matrix) const
{
    const Eigen::Index outerCount = matrix.outerSize() + 1;
    const Eigen::Index nonZeroCount = matrix.nonZeros();
    return initialized_
        && dimension_ == matrix.rows()
        && outerIndices_.size() == static_cast<std::size_t>(outerCount)
        && innerIndices_.size() == static_cast<std::size_t>(nonZeroCount)
        && arraysEqual(outerIndices_.data(), matrix.outerIndexPtr(), outerCount)
        && arraysEqual(innerIndices_.data(), matrix.innerIndexPtr(), nonZeroCount);
}

void PardisoSolver::configureThreads(int threadCount)
{
    if (threadCount < 0) {
        throw std::invalid_argument("PARDISO thread count must be non-negative");
    }
    activeThreadCount_ = threadCount > 0 ? threadCount : defaultThreadCount_;
}

void PardisoSolver::runPhase(
    int phase,
    const Eigen::SparseMatrix<double>& matrix,
    double* rightHandSide,
    double* result)
{
    const std::size_t maximumParallelism = static_cast<std::size_t>(
        (std::max)(1, activeThreadCount_));
    tbb::global_control threadLimit(
        tbb::global_control::max_allowed_parallelism,
        maximumParallelism);
    const auto begin = std::chrono::steady_clock::now();
    const int maxFactorizations = 1;
    const int matrixNumber = 1;
    const int matrixType = realSymmetricPositiveDefinite;
    const int rightHandSideCount = 1;
    const int messageLevel = 0;
    int error = 0;

    // A lower-triangular CSC matrix is the same memory layout as the upper
    // triangle of its symmetric transpose in CSR, which is what PARDISO uses.
    pardiso(
        internalData_.data(),
        &maxFactorizations,
        &matrixNumber,
        &matrixType,
        &phase,
        &dimension_,
        matrix.valuePtr(),
        matrix.outerIndexPtr(),
        matrix.innerIndexPtr(),
        permutation_.data(),
        &rightHandSideCount,
        parameters_.data(),
        &messageLevel,
        rightHandSide,
        result,
        &error);
    const double elapsedMilliseconds = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
    if (phase == analysisPhase) {
        analysisMilliseconds_ += elapsedMilliseconds;
    }
    else if (phase == factorizationPhase) {
        factorizationMilliseconds_ += elapsedMilliseconds;
    }
    else if (phase == solvePhase) {
        solveMilliseconds_ += elapsedMilliseconds;
    }
    if (error != 0) {
        throw std::runtime_error(
            "oneMKL PARDISO phase " + std::to_string(phase)
            + " failed with error " + std::to_string(error));
    }
}

void PardisoSolver::solve(
    const Eigen::SparseMatrix<double>& lowerTriangle,
    const Eigen::VectorXd& rightHandSide,
    Eigen::VectorXd& result,
    int threadCount)
{
    if (lowerTriangle.rows() != lowerTriangle.cols()) {
        throw std::invalid_argument("PARDISO requires a square matrix");
    }
    if (!lowerTriangle.isCompressed()) {
        throw std::invalid_argument("PARDISO requires compressed sparse input");
    }
    if (rightHandSide.size() != lowerTriangle.rows()) {
        throw std::invalid_argument("PARDISO right-hand side has an incompatible size");
    }

    configureThreads(threadCount);
    if (!hasSamePattern(lowerTriangle) || forceReordering_) {
        release();
        initialize();
        dimension_ = static_cast<int>(lowerTriangle.rows());
        permutation_.resize(dimension_);
        const bool reusePermutation = permutationReady_
            && permutationDimension_ == dimension_
            && !forceReordering_;
        forceReordering_ = false;
        updateFactorReference_ = !reusePermutation;
        parameters_[4] = reusePermutation ? 1 : 2;
        initialized_ = true;
        runPhase(analysisPhase, lowerTriangle, nullptr, nullptr);
        parameters_[4] = 0;
        if (!reusePermutation) {
            permutationReady_ = true;
            permutationDimension_ = dimension_;
        }
        outerIndices_.assign(
            lowerTriangle.outerIndexPtr(),
            lowerTriangle.outerIndexPtr() + lowerTriangle.outerSize() + 1);
        innerIndices_.assign(
            lowerTriangle.innerIndexPtr(),
            lowerTriangle.innerIndexPtr() + lowerTriangle.nonZeros());
        ++symbolicAnalysisCount_;
    }

    runPhase(factorizationPhase, lowerTriangle, nullptr, nullptr);
    factorNonZeros_ = parameters_[17];
    if (updateFactorReference_) {
        referenceFactorNonZeros_ = factorNonZeros_;
        updateFactorReference_ = false;
    }
    else if (referenceFactorNonZeros_ > 0
        && factorNonZeros_ > static_cast<int>(1.2 * referenceFactorNonZeros_)) {
        forceReordering_ = true;
    }
    ++numericFactorizationCount_;

    rightHandSideBuffer_ = rightHandSide;
    result.resize(rightHandSide.size());
    runPhase(
        solvePhase,
        lowerTriangle,
        rightHandSideBuffer_.data(),
        result.data());
    if (!result.allFinite()) {
        throw std::runtime_error("oneMKL PARDISO produced a non-finite result");
    }
}
