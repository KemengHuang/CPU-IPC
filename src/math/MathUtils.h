#pragma once
#ifndef FEM_MATH_H
#define FEM_MATH_H

#include "fem/SVD.h"
#include <vector>

// Sum of squared norms over a collection of 3-vectors.
double vector_squareNorm(const std::vector<Eigen::Vector3d>& vecs);

// ------------------------------------------------------------------
// Dynamic-size vectorization helpers (backward compatible)
// Stack the entries of a matrix column by column into a vector.
// ------------------------------------------------------------------
Eigen::MatrixXd vec_double(const Eigen::MatrixXd& F);
Eigen::MatrixXf vec_float(const Eigen::MatrixXf& F);

// ------------------------------------------------------------------
// Fixed-size vectorization helpers for hot paths
// These avoid heap allocations for the common 3x3, 3x2 and 3x4 shapes.
// ------------------------------------------------------------------
Eigen::Matrix<double, 9, 1> vec9(const Eigen::Matrix3d& F);
Eigen::Matrix<double, 6, 1> vec6(const Eigen::Matrix<double, 3, 2>& F);
Eigen::Matrix<double, 12, 1> vec12(const Eigen::Matrix<double, 3, 4>& F);

// ------------------------------------------------------------------
// Cubic equation solvers
// All solve a*x^3 + b*x^2 + c*x + d = 0.
// ------------------------------------------------------------------

// Newton-like iterative solver used by the Stable Neo-Hookean eigen-
// projection.  Returns three roots (possibly coincident).
std::vector<double> NewtonSolverForCubicEquation(const double& a, const double& b, const double& c, const double& d);

// Variant cubic solver used for degenerate / special-case handling.
std::vector<double> NewtonSolverForCubicEquation_snk(const double& a, const double& b, const double& c, const double& d);

// Closed-form cubic solver (Vieta / Cardano style).
std::vector<double> __SolverForCubicEquation(const double& a, const double& b, const double& c, const double& d, double EPS = 1e-6);

#endif // !FEM_MATH_H
