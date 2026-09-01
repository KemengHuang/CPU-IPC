#pragma once

#include <Eigen/Core>

#include <vector>

namespace ElasticityMath {

Eigen::MatrixXd vectorizeColumnMajor(const Eigen::MatrixXd& matrix);

std::vector<double> solveStableNeoHookeanCubic(
    double a,
    double b,
    double c,
    double d);

std::vector<double> solveCubicRealRoots(
    double a,
    double b,
    double c,
    double d,
    double tolerance = 1e-6);

} // namespace ElasticityMath
