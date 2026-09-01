#include "SuiteSparseLDLSolver.h"

#include <amd.h>
#include <ldl.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace {

template <typename Scalar>
bool arraysEqual(const Scalar* lhs, const Scalar* rhs, Eigen::Index count)
{
    return count == 0
        || std::memcmp(lhs, rhs, static_cast<std::size_t>(count) * sizeof(Scalar)) == 0;
}

} // namespace

bool SuiteSparseLDLSolver::hasSamePattern(
    const Eigen::SparseMatrix<double>& lowerTriangle) const
{
    const Eigen::Index outerCount = lowerTriangle.outerSize() + 1;
    const Eigen::Index nonZeroCount = lowerTriangle.nonZeros();
    return patternReady_
        && outerIndices_.size() == static_cast<std::size_t>(outerCount)
        && innerIndices_.size() == static_cast<std::size_t>(nonZeroCount)
        && arraysEqual(outerIndices_.data(), lowerTriangle.outerIndexPtr(), outerCount)
        && arraysEqual(innerIndices_.data(), lowerTriangle.innerIndexPtr(), nonZeroCount);
}

void SuiteSparseLDLSolver::analyzePattern(
    const Eigen::SparseMatrix<double>& lowerTriangle)
{
    const int dimension = static_cast<int>(lowerTriangle.rows());
    permutation_.resize(dimension);
    inversePermutation_.resize(dimension);

    if (dimension % 3 == 0) {
        const int blockDimension = dimension / 3;
        std::vector<Eigen::Triplet<double>> blockTriplets;
        blockTriplets.reserve(lowerTriangle.nonZeros());
        for (int column = 0; column < lowerTriangle.outerSize(); ++column) {
            for (Eigen::SparseMatrix<double>::InnerIterator entry(
                    lowerTriangle, column); entry; ++entry) {
                blockTriplets.emplace_back(entry.row() / 3, entry.col() / 3, 1.0);
            }
        }

        Eigen::SparseMatrix<double> blockGraph(blockDimension, blockDimension);
        blockGraph.setFromTriplets(
            blockTriplets.begin(), blockTriplets.end(),
            [](double, double) { return 1.0; });
        std::vector<int> blockPermutation(blockDimension);
        if (amd_order(
                blockDimension,
                blockGraph.outerIndexPtr(),
                blockGraph.innerIndexPtr(),
                blockPermutation.data(), nullptr, nullptr) < AMD_OK) {
            throw std::runtime_error("SuiteSparse block AMD ordering failed");
        }
        for (int block = 0; block < blockDimension; ++block) {
            for (int component = 0; component < 3; ++component) {
                permutation_[3 * block + component]
                    = 3 * blockPermutation[block] + component;
            }
        }
    }
    else if (amd_order(
                 dimension,
                 lowerTriangle.outerIndexPtr(),
                 lowerTriangle.innerIndexPtr(),
                 permutation_.data(), nullptr, nullptr) < AMD_OK) {
        throw std::runtime_error("SuiteSparse AMD ordering failed");
    }
    for (int permuted = 0; permuted < dimension; ++permuted) {
        inversePermutation_[permutation_[permuted]] = permuted;
    }

    std::vector<Eigen::Triplet<double>> permutedTriplets;
    permutedTriplets.reserve(lowerTriangle.nonZeros());
    for (int column = 0; column < lowerTriangle.outerSize(); ++column) {
        for (Eigen::SparseMatrix<double>::InnerIterator entry(
                lowerTriangle, column); entry; ++entry) {
            if (entry.row() < entry.col()) {
                throw std::invalid_argument(
                    "SuiteSparse LDL input contains entries above the diagonal");
            }
            int row = inversePermutation_[entry.row()];
            int permutedColumn = inversePermutation_[entry.col()];
            if (row > permutedColumn) {
                (std::swap)(row, permutedColumn);
            }
            permutedTriplets.emplace_back(row, permutedColumn, entry.value());
        }
    }
    matrix_.resize(dimension, dimension);
    matrix_.setFromTriplets(permutedTriplets.begin(), permutedTriplets.end());
    matrix_.makeCompressed();
    if (matrix_.nonZeros() != lowerTriangle.nonZeros()) {
        throw std::logic_error("SuiteSparse LDL permutation changed the sparse entry count");
    }
    if (!ldl_valid_matrix(
            dimension,
            matrix_.outerIndexPtr(),
            matrix_.innerIndexPtr())) {
        throw std::runtime_error("SuiteSparse LDL produced an invalid permuted matrix");
    }

    valueIndices_.resize(lowerTriangle.nonZeros());
    Eigen::Index source = 0;
    for (int column = 0; column < lowerTriangle.outerSize(); ++column) {
        for (Eigen::SparseMatrix<double>::InnerIterator entry(
                lowerTriangle, column); entry; ++entry, ++source) {
            int row = inversePermutation_[entry.row()];
            int permutedColumn = inversePermutation_[entry.col()];
            if (row > permutedColumn) {
                (std::swap)(row, permutedColumn);
            }

            const int begin = matrix_.outerIndexPtr()[permutedColumn];
            const int end = matrix_.outerIndexPtr()[permutedColumn + 1];
            const int* position = std::lower_bound(
                matrix_.innerIndexPtr() + begin,
                matrix_.innerIndexPtr() + end,
                row);
            if (position == matrix_.innerIndexPtr() + end || *position != row) {
                throw std::logic_error("permuted LDL entry is missing from the sparse matrix");
            }
            valueIndices_[source] = position - matrix_.innerIndexPtr();
        }
    }

    columnPointers_.resize(dimension + 1);
    eliminationTree_.resize(dimension);
    columnNonZeros_.resize(dimension);
    flags_.resize(dimension);
    ldl_symbolic(
        dimension,
        matrix_.outerIndexPtr(),
        matrix_.innerIndexPtr(),
        columnPointers_.data(),
        eliminationTree_.data(),
        columnNonZeros_.data(),
        flags_.data(),
        nullptr,
        nullptr);

    const int factorNonZeros = columnPointers_.back();
    if (factorNonZeros < 0) {
        throw std::runtime_error("SuiteSparse LDL symbolic analysis returned an invalid size");
    }
    rowIndices_.resize(factorNonZeros);
    factorValues_.resize(factorNonZeros);
    diagonal_.resize(dimension);
    numericWorkspace_.resize(dimension);
    patternWorkspace_.resize(dimension);

    outerIndices_.assign(
        lowerTriangle.outerIndexPtr(),
        lowerTriangle.outerIndexPtr() + lowerTriangle.outerSize() + 1);
    innerIndices_.assign(
        lowerTriangle.innerIndexPtr(),
        lowerTriangle.innerIndexPtr() + lowerTriangle.nonZeros());
    patternReady_ = true;
    ++symbolicAnalysisCount_;
}

void SuiteSparseLDLSolver::updateValues(
    const Eigen::SparseMatrix<double>& lowerTriangle)
{
    if (valueIndices_.size() != static_cast<std::size_t>(lowerTriangle.nonZeros())) {
        throw std::logic_error("SuiteSparse LDL value map has an incompatible size");
    }

    const double* sourceValues = lowerTriangle.valuePtr();
    double* destinationValues = matrix_.valuePtr();
    for (Eigen::Index source = 0; source < lowerTriangle.nonZeros(); ++source) {
        destinationValues[valueIndices_[source]] = sourceValues[source];
    }
}

void SuiteSparseLDLSolver::factorize()
{
    const int dimension = static_cast<int>(matrix_.rows());
    if (ldl_numeric(
            dimension,
            matrix_.outerIndexPtr(),
            matrix_.innerIndexPtr(),
            matrix_.valuePtr(),
            columnPointers_.data(),
            eliminationTree_.data(),
            columnNonZeros_.data(),
            rowIndices_.data(),
            factorValues_.data(),
            diagonal_.data(),
            numericWorkspace_.data(),
            patternWorkspace_.data(),
            flags_.data(),
            nullptr,
            nullptr) != dimension) {
        throw std::runtime_error("SuiteSparse LDL factorization failed");
    }

    for (double diagonalValue : diagonal_) {
        if (!std::isfinite(diagonalValue) || diagonalValue <= 0.0) {
            throw std::runtime_error("SuiteSparse LDL factorization is not positive definite");
        }
    }
    ++numericFactorizationCount_;
}

void SuiteSparseLDLSolver::solveFactorized(
    const Eigen::VectorXd& rightHandSide,
    Eigen::VectorXd& result)
{
    const int dimension = static_cast<int>(matrix_.rows());
    if (rightHandSide.size() != dimension) {
        throw std::invalid_argument("SuiteSparse LDL right-hand side has an incompatible size");
    }

    for (int permuted = 0; permuted < dimension; ++permuted) {
        numericWorkspace_[permuted] = rightHandSide[permutation_[permuted]];
    }
    ldl_lsolve(
        dimension,
        numericWorkspace_.data(),
        columnPointers_.data(),
        rowIndices_.data(),
        factorValues_.data());
    ldl_dsolve(dimension, numericWorkspace_.data(), diagonal_.data());
    ldl_ltsolve(
        dimension,
        numericWorkspace_.data(),
        columnPointers_.data(),
        rowIndices_.data(),
        factorValues_.data());

    result.resize(dimension);
    for (int permuted = 0; permuted < dimension; ++permuted) {
        result[permutation_[permuted]] = numericWorkspace_[permuted];
    }
    if (!result.allFinite()) {
        throw std::runtime_error("SuiteSparse LDL solve produced a non-finite result");
    }
}

void SuiteSparseLDLSolver::solve(
    const Eigen::SparseMatrix<double>& lowerTriangle,
    const Eigen::VectorXd& rightHandSide,
    Eigen::VectorXd& result)
{
    if (lowerTriangle.rows() != lowerTriangle.cols()) {
        throw std::invalid_argument("SuiteSparse LDL requires a square matrix");
    }
    if (!lowerTriangle.isCompressed()) {
        throw std::invalid_argument("SuiteSparse LDL requires compressed sparse input");
    }

    if (!hasSamePattern(lowerTriangle)) {
        analyzePattern(lowerTriangle);
    }
    else {
        updateValues(lowerTriangle);
    }
    factorize();
    solveFactorized(rightHandSide, result);
}
