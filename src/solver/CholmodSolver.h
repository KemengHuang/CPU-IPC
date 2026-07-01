//
// Created by lamws on 2021/8/3.
//

#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <suitesparse/cholmod.h>

class CholmodSolver {
public:
    CholmodSolver();

    ~CholmodSolver();

    void set_pattern(const Eigen::SparseMatrix<double>& mtr);

    void set_pattern(const std::vector<Eigen::Triplet<double>>& input, int nrow, int ncol);

    void solve(Eigen::VectorXd& rhs,
        Eigen::VectorXd& result);

    void preFactorize(const Eigen::SparseMatrix<double>& mtr);
    void solve_with_preFactorize(Eigen::VectorXd& rhs, Eigen::VectorXd& result);
protected:
    int numRows;
    Eigen::VectorXi ia, ja;
    std::vector<std::map<int, int>> IJ2aI;
    Eigen::VectorXd a;

    cholmod_triplet* triplets;
    cholmod_common cm;
    cholmod_sparse* A;
    cholmod_factor* L;
    cholmod_dense* b, * solution;
    cholmod_dense* x_cd, * y_cd; // for multiply

    void* Ai, * Ap, * Ax, * bx, * solutionx, * x_cdx, * y_cdx;
};
