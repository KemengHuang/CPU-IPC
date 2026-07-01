#pragma once
#ifndef FEM_MATH_H
#define FEM_MATH_H

#include "SVD.h"
#include<vector>

double vector_squareNorm(std::vector<Eigen::Vector3d> vecs);
Eigen::MatrixXd vec_double(const Eigen::MatrixXd& F);
Eigen::MatrixXf vec_float(const Eigen::MatrixXf& F);
std::vector<double> NewtonSolverForCubicEquation(const double& a, const double& b, const double& c, const double& d);
std::vector<double> NewtonSolverForCubicEquation_snk(const double& a, const double& b, const double& c, const double& d);
bool segmentHitTest(const Eigen::Vector3d& point0, const Eigen::Vector3d& point1, const Eigen::Vector3d& point2, const Eigen::Vector3d& seg_point0, const Eigen::Vector3d& seg_point1);
std::vector<double> __SolverForCubicEquation(const double& a, const double& b, const double& c, const double& d, double EPS = 1e-6);
#endif // !FEM_MATH_H
