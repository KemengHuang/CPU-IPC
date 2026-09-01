#pragma once

#include <Eigen/Core>

#include <vector>

struct HingeBendingInfo {
    Eigen::Vector4i vertices;
    double restAngle = 0.0;
    double geometricWeight = 0.0;
};

HingeBendingInfo makeHingeBendingInfo(
    int edgeStart,
    int edgeEnd,
    int firstOpposite,
    int secondOpposite,
    const std::vector<Eigen::Vector3d>& restPositions);

double hingeBendingEnergy(
    const HingeBendingInfo& hinge,
    const std::vector<Eigen::Vector3d>& positions,
    double plateRigidity);

void hingeBendingGradientAndHessian(
    const HingeBendingInfo& hinge,
    const std::vector<Eigen::Vector3d>& positions,
    double plateRigidity,
    Eigen::Matrix<double, 12, 1>& gradient,
    Eigen::Matrix<double, 12, 12>& hessian);
