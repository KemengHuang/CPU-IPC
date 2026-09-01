#include "TetInversionGuard.h"

#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

namespace {

double determinant(
    const Eigen::Vector3d& first,
    const Eigen::Vector3d& second,
    const Eigen::Vector3d& third)
{
    return first.dot(second.cross(third));
}

double smallestPositiveQuadraticRoot(
    double a,
    double b,
    double c,
    double tolerance)
{
    if (std::abs(a) <= tolerance) {
        if (std::abs(b) <= tolerance) {
            return -1.0;
        }
        const double root = -c / b;
        return root > 0.0 ? root : -1.0;
    }

    const double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0) {
        return -1.0;
    }
    const double squareRoot = std::sqrt(discriminant);
    const double first = (-b - squareRoot) / (2.0 * a);
    const double second = (-b + squareRoot) / (2.0 * a);
    double result = std::numeric_limits<double>::infinity();
    if (first > 0.0) {
        result = first;
    }
    if (second > 0.0) {
        result = (std::min)(result, second);
    }
    return std::isfinite(result) ? result : -1.0;
}

double smallestPositiveCubicRoot(
    double a,
    double b,
    double c,
    double d,
    double tolerance)
{
    if (std::abs(a) <= tolerance) {
        return smallestPositiveQuadraticRoot(b, c, d, tolerance);
    }

    using Complex = std::complex<double>;
    const Complex delta0(b * b - 3.0 * a * c, 0.0);
    const Complex delta1(
        2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d,
        0.0);
    const Complex discriminant = delta1 * delta1
        - 4.0 * delta0 * delta0 * delta0;
    Complex cubeRoot = std::pow(
        (delta1 + std::sqrt(discriminant)) / 2.0,
        1.0 / 3.0);
    if (std::abs(cubeRoot) <= tolerance) {
        cubeRoot = std::pow(
            (delta1 - std::sqrt(discriminant)) / 2.0,
            1.0 / 3.0);
    }
    if (std::abs(cubeRoot) <= tolerance) {
        const double repeatedRoot = -b / (3.0 * a);
        return repeatedRoot > 0.0 ? repeatedRoot : -1.0;
    }

    const Complex omega(-0.5, std::sqrt(3.0) / 2.0);
    const Complex roots[3] = {
        -(b + cubeRoot + delta0 / cubeRoot) / (3.0 * a),
        -(b + omega * cubeRoot + delta0 / (omega * cubeRoot)) / (3.0 * a),
        -(b + std::conj(omega) * cubeRoot
            + delta0 / (std::conj(omega) * cubeRoot)) / (3.0 * a)
    };

    double result = std::numeric_limits<double>::infinity();
    for (const Complex& root : roots) {
        if (std::abs(root.imag()) <= tolerance && root.real() > 0.0) {
            result = (std::min)(result, root.real());
        }
    }
    return std::isfinite(result) ? result : -1.0;
}

double tetVolumeLimit(
    const mesh3D& mesh,
    const Eigen::Vector4i& tet,
    const std::vector<Eigen::Vector3d>& searchDirection,
    double volumeRatio)
{
    const Eigen::Vector3d first = mesh.vertexes[tet[1]] - mesh.vertexes[tet[0]];
    const Eigen::Vector3d second = mesh.vertexes[tet[2]] - mesh.vertexes[tet[0]];
    const Eigen::Vector3d third = mesh.vertexes[tet[3]] - mesh.vertexes[tet[0]];
    const Eigen::Vector3d firstMotion =
        -searchDirection[tet[1]] + searchDirection[tet[0]];
    const Eigen::Vector3d secondMotion =
        -searchDirection[tet[2]] + searchDirection[tet[0]];
    const Eigen::Vector3d thirdMotion =
        -searchDirection[tet[3]] + searchDirection[tet[0]];

    const double constant = (1.0 - volumeRatio)
        * determinant(first, second, third);
    const double linear =
        determinant(firstMotion, second, third)
        + determinant(first, secondMotion, third)
        + determinant(first, second, thirdMotion);
    const double quadratic =
        determinant(firstMotion, secondMotion, third)
        + determinant(firstMotion, second, thirdMotion)
        + determinant(first, secondMotion, thirdMotion);
    const double cubic = determinant(
        firstMotion, secondMotion, thirdMotion);
    return smallestPositiveCubicRoot(
        cubic, quadratic, linear, constant, 1e-6);
}

} // namespace

void limitStepToPreventTetInversion(
    const mesh3D& mesh,
    const std::vector<Eigen::Vector3d>& searchDirection,
    double volumeRatio,
    double& stepSize)
{
    if (mesh.tetrahedras.empty()) {
        return;
    }
    if (searchDirection.size() != mesh.vertexes.size()
        || volumeRatio <= 0.0 || volumeRatio >= 1.0) {
        throw std::invalid_argument("invalid tetrahedron inversion-guard input");
    }

    Eigen::VectorXd limits(mesh.tetrahedras.size());
    tbb::parallel_for(0, static_cast<int>(mesh.tetrahedras.size()), 1,
        [&](int tet) {
            const double root = tetVolumeLimit(
                mesh, mesh.tetrahedras[tet], searchDirection, volumeRatio);
            limits[tet] = root > 0.0
                ? root
                : std::numeric_limits<double>::infinity();
        });
    stepSize = (std::min)(stepSize, limits.minCoeff());
}
