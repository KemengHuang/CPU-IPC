#include "HingeBending.h"

#include "Elasticity.h"

#include <cmath>
#include <stdexcept>

namespace {

void validateVertexIndex(
    int vertex,
    const std::vector<Eigen::Vector3d>& positions)
{
    if (vertex < 0 || vertex >= static_cast<int>(positions.size())) {
        throw std::out_of_range("hinge bending vertex index is out of range");
    }
}

double angleForHinge(
    const HingeBendingInfo& hinge,
    const std::vector<Eigen::Vector3d>& positions,
    Eigen::Matrix<double, 1, 12>* gradient,
    Eigen::Matrix<double, 12, 12>* hessian)
{
    for (int corner = 0; corner < 4; ++corner) {
        validateVertexIndex(hinge.vertices[corner], positions);
    }
    return edgeTheta(
        positions[hinge.vertices[0]],
        positions[hinge.vertices[1]],
        positions[hinge.vertices[2]],
        positions[hinge.vertices[3]],
        gradient,
        hessian);
}

double angleDifference(double angle, double restAngle)
{
    constexpr double twoPi = 6.28318530717958647692;
    return std::remainder(angle - restAngle, twoPi);
}

} // namespace

HingeBendingInfo makeHingeBendingInfo(
    int edgeStart,
    int edgeEnd,
    int firstOpposite,
    int secondOpposite,
    const std::vector<Eigen::Vector3d>& restPositions)
{
    HingeBendingInfo result;
    result.vertices = Eigen::Vector4i(
        edgeStart, edgeEnd, firstOpposite, secondOpposite);
    for (int corner = 0; corner < 4; ++corner) {
        validateVertexIndex(result.vertices[corner], restPositions);
    }

    const Eigen::Vector3d edge =
        restPositions[edgeEnd] - restPositions[edgeStart];
    const double edgeLength = edge.norm();
    const double twiceFirstArea = edge.cross(
        restPositions[firstOpposite] - restPositions[edgeStart]).norm();
    const double twiceSecondArea = edge.cross(
        restPositions[secondOpposite] - restPositions[edgeStart]).norm();
    const double twiceAreaSum = twiceFirstArea + twiceSecondArea;
    if (!std::isfinite(edgeLength) || !std::isfinite(twiceAreaSum)
        || edgeLength <= 0.0 || twiceAreaSum <= 0.0) {
        throw std::runtime_error("degenerate rest hinge geometry");
    }

    // h0 + h1 = 2(A0 + A1) / l0, hence l0/(h0+h1).
    result.geometricWeight = edgeLength * edgeLength / twiceAreaSum;
    result.restAngle = angleForHinge(result, restPositions, nullptr, nullptr);
    return result;
}

double hingeBendingEnergy(
    const HingeBendingInfo& hinge,
    const std::vector<Eigen::Vector3d>& positions,
    double plateRigidity)
{
    const double angle = angleForHinge(hinge, positions, nullptr, nullptr);
    const double difference = angleDifference(angle, hinge.restAngle);
    return 0.5 * plateRigidity * hinge.geometricWeight
        * difference * difference;
}

void hingeBendingGradientAndHessian(
    const HingeBendingInfo& hinge,
    const std::vector<Eigen::Vector3d>& positions,
    double plateRigidity,
    Eigen::Matrix<double, 12, 1>& gradient,
    Eigen::Matrix<double, 12, 12>& hessian)
{
    Eigen::Matrix<double, 1, 12> angleGradient;
    Eigen::Matrix<double, 12, 12> angleHessian;
    const double angle = angleForHinge(
        hinge, positions, &angleGradient, &angleHessian);
    const double difference = angleDifference(angle, hinge.restAngle);
    const double scale = plateRigidity * hinge.geometricWeight;
    gradient.noalias() = scale * difference * angleGradient.transpose();
    hessian.noalias() = scale * (
        angleGradient.transpose() * angleGradient
        + difference * angleHessian);
}
