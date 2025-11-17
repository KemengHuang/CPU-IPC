#pragma once
#ifndef FEM_MATH_H
#define FEM_MATH_H

#include "SVD.h"
#include<vector>

double vector_squareNorm(std::vector<Vector3d> vecs);
MatrixXd vec_double(MatrixXd F);
MatrixXf vec_float(MatrixXf F);
std::vector<double> NewtonSolverForCubicEquation(const double& a, const double& b, const double& c, const double& d);
std::vector<double> NewtonSolverForCubicEquation_snk(const double& a, const double& b, const double& c, const double& d);
bool segmentHitTest(const Vector3d& point0, const Vector3d& point1, const Vector3d& point2, const Vector3d& seg_point0, const Vector3d& seg_point1);
std::vector<double> __SolverForCubicEquation(const double& a, const double& b, const double& c, const double& d, double EPS = 1e-6);
#endif // !FEM_MATH_H

