#include "CholmodSolver.h"

#include <cstring>
#include <stdexcept>

namespace {

template <typename Scalar>
bool arraysEqual(const Scalar* lhs, const Scalar* rhs, Eigen::Index count)
{
    return count == 0 || std::memcmp(lhs, rhs, static_cast<size_t>(count) * sizeof(Scalar)) == 0;
}

void configureSymbolicOrdering(cholmod_common& common)
{
#ifdef CIPC_CHOLMOD_USE_METIS
    // Keep CHOLMOD's cost-aware automatic policy: evaluate AMD first, then use
    // METIS when AMD predicts excessive fill or work. Explicitly select METIS
    // (rather than NESDIS) as the nested-dissection fallback.
    common.nmethods = 0;
    common.default_nesdis = 0;
#else
    // Make the build option a real A/B switch instead of leaving METIS
    // reachable through CHOLMOD's default nmethods == 0 policy.
    common.nmethods = 1;
    common.method[0].ordering = CHOLMOD_AMD;
#endif
    common.postorder = 1;
}

} // namespace

CholmodSolver::CholmodSolver()
{
    std::memset(&common_, 0, sizeof(common_));
    std::memset(&matrix_, 0, sizeof(matrix_));
    if (!cholmod_start(&common_)) {
        throw std::runtime_error("CHOLMOD initialization failed");
    }
    configureSymbolicOrdering(common_);
}

CholmodSolver::~CholmodSolver()
{
    if (factor_ != nullptr) {
        cholmod_free_factor(&factor_, &common_);
    }
    cholmod_finish(&common_);
}

bool CholmodSolver::hasSamePattern(const Eigen::SparseMatrix<double>& matrix) const
{
    if (!matrixReady_ || numRows_ != matrix.rows() || numCols_ != matrix.cols()) {
        return false;
    }

    const Eigen::Index outerCount = matrix.outerSize() + 1;
    const Eigen::Index nonZeroCount = matrix.nonZeros();
    if (outerIndices_.size() != outerCount || innerIndices_.size() != nonZeroCount) {
        return false;
    }

    return arraysEqual(outerIndices_.data(), matrix.outerIndexPtr(), outerCount)
        && arraysEqual(innerIndices_.data(), matrix.innerIndexPtr(), nonZeroCount);
}

void CholmodSolver::updateMatrixView()
{
    matrix_.nrow = static_cast<size_t>(numRows_);
    matrix_.ncol = static_cast<size_t>(numCols_);
    matrix_.nzmax = static_cast<size_t>(values_.size());
    matrix_.p = outerIndices_.data();
    matrix_.i = innerIndices_.data();
    matrix_.x = values_.data();
    matrix_.z = nullptr;
    matrix_.nz = nullptr;
    matrix_.stype = -1;
    matrix_.itype = CHOLMOD_INT;
    matrix_.xtype = CHOLMOD_REAL;
    matrix_.dtype = CHOLMOD_DOUBLE;
    matrix_.sorted = true;
    matrix_.packed = true;
}

void CholmodSolver::set_pattern(const Eigen::SparseMatrix<double>& input)
{
    if (input.rows() != input.cols()) {
        throw std::invalid_argument("CHOLMOD requires a square matrix");
    }

    if (!input.isCompressed()) {
        Eigen::SparseMatrix<double> compressed = input;
        compressed.makeCompressed();
        set_pattern(compressed);
        return;
    }

    const bool patternChanged = !hasSamePattern(input);
    if (patternChanged && factor_ != nullptr) {
        cholmod_free_factor(&factor_, &common_);
    }

    numRows_ = static_cast<int>(input.rows());
    numCols_ = static_cast<int>(input.cols());

    outerIndices_.resize(input.outerSize() + 1);
    std::memcpy(
        outerIndices_.data(),
        input.outerIndexPtr(),
        static_cast<size_t>(outerIndices_.size()) * sizeof(outerIndices_[0]));

    innerIndices_.resize(input.nonZeros());
    if (input.nonZeros() > 0) {
        std::memcpy(
            innerIndices_.data(),
            input.innerIndexPtr(),
            static_cast<size_t>(innerIndices_.size()) * sizeof(innerIndices_[0]));
    }

    values_.resize(input.nonZeros());
    if (input.nonZeros() > 0) {
        std::memcpy(
            values_.data(),
            input.valuePtr(),
            static_cast<size_t>(values_.size()) * sizeof(values_[0]));
    }

    updateMatrixView();
    matrixReady_ = true;
}

void CholmodSolver::factorize()
{
    if (!matrixReady_) {
        throw std::logic_error("CHOLMOD matrix has not been initialized");
    }

    if (factor_ == nullptr) {
        factor_ = cholmod_analyze(&matrix_, &common_);
        if (factor_ == nullptr) {
            throw std::runtime_error("CHOLMOD symbolic analysis failed");
        }
        ++symbolicAnalysisCount_;
    }

    if (!cholmod_factorize(&matrix_, factor_, &common_)) {
        throw std::runtime_error("CHOLMOD numeric factorization failed");
    }
    ++numericFactorizationCount_;
}

void CholmodSolver::solveFactorized(const Eigen::VectorXd& rhs, Eigen::VectorXd& result)
{
    if (factor_ == nullptr) {
        throw std::logic_error("CHOLMOD matrix has not been factorized");
    }
    if (rhs.size() != numRows_) {
        throw std::invalid_argument("CHOLMOD right-hand side has an incompatible size");
    }

    cholmod_dense rhsView;
    std::memset(&rhsView, 0, sizeof(rhsView));
    rhsView.nrow = static_cast<size_t>(rhs.size());
    rhsView.ncol = 1;
    rhsView.nzmax = static_cast<size_t>(rhs.size());
    rhsView.d = static_cast<size_t>(rhs.size());
    rhsView.x = const_cast<double*>(rhs.data());
    rhsView.xtype = CHOLMOD_REAL;
    rhsView.dtype = CHOLMOD_DOUBLE;

    cholmod_dense* solution = cholmod_solve(CHOLMOD_A, factor_, &rhsView, &common_);
    if (solution == nullptr) {
        throw std::runtime_error("CHOLMOD solve failed");
    }

    result.resize(rhs.size());
    std::memcpy(
        result.data(),
        solution->x,
        static_cast<size_t>(result.size()) * sizeof(result[0]));
    cholmod_free_dense(&solution, &common_);
}

void CholmodSolver::solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& result)
{
    factorize();
    solveFactorized(rhs, result);
}

void CholmodSolver::preFactorize(const Eigen::SparseMatrix<double>& matrix)
{
    set_pattern(matrix);
    factorize();
}

void CholmodSolver::solve_with_preFactorize(const Eigen::VectorXd& rhs, Eigen::VectorXd& result)
{
    solveFactorized(rhs, result);
}
