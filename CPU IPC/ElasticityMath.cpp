#include "ElasticityMath.h"

#include <cmath>

namespace {

double evaluateCubic(double x, double a, double b, double c, double d)
{
    return a * x * x * x + b * x * x + c * x + d;
}

double evaluateCubicDerivative(double x, double a, double b, double c)
{
    return 3 * a * x * x + 2 * b * x + c;
}

} // namespace

namespace ElasticityMath {

Eigen::MatrixXd vectorizeColumnMajor(const Eigen::MatrixXd& matrix)
{
    Eigen::MatrixXd result(matrix.size(), 1);
    for (int column = 0; column < matrix.cols(); ++column) {
        for (int row = 0; row < matrix.rows(); ++row) {
            result(column * matrix.rows() + row, 0) = matrix(row, column);
        }
    }
    return result;
}

std::vector<double> solveStableNeoHookeanCubic(
    double a,
    double b,
    double c,
    double d)
{
    constexpr double tolerance = 1e-6;
    double step = 0.0;
    std::vector<double> roots(3);
    const double specialPoint = -b / a / 3;
    double stationaryPoints[2];
    int rootCount = 1;
    const double discriminant = 4 * b * b - 12 * a * c;
    int resultIndex = 0;
    if (discriminant > 0) {
        stationaryPoints[0] = (std::sqrt(discriminant) - 2 * b) / (6 * a);
        stationaryPoints[1] = (-std::sqrt(discriminant) - 2 * b) / (6 * a);
        double firstValue = evaluateCubic(stationaryPoints[0], a, b, c, d);
        double secondValue = evaluateCubic(stationaryPoints[1], a, b, c, d);
        if (firstValue >= 0) {
            firstValue = 0;
            roots[1] = stationaryPoints[0];
            roots[2] = stationaryPoints[0];
        }
        if (secondValue <= 0) {
            secondValue = 0;
            roots[1] = stationaryPoints[1];
            roots[2] = stationaryPoints[1];
        }
        step = stationaryPoints[0] - stationaryPoints[1];
        if (firstValue * secondValue < 0) {
            rootCount = 3;
        }
        else if (secondValue <= 0) {
            step = -step;
        }
    }
    else {
        roots[0] = specialPoint;
        roots[1] = specialPoint;
        roots[2] = specialPoint;
        return roots;
    }

    double start = specialPoint - step;
    double previous = start;
    for (int root = 0; root < rootCount; ++root) {
        double current = 0.0;
        int iteration = 0;
        do {
            if (iteration > 0) {
                previous = current;
            }
            current = previous
                - evaluateCubic(previous, a, b, c, d)
                    / evaluateCubicDerivative(previous, a, b, c);
            ++iteration;
        } while (std::abs(current - previous) > tolerance * tolerance
            && iteration < 100000);
        roots[resultIndex++] = current;
        start += step;
        previous = start;
    }
    return roots;
}

std::vector<double> solveCubicRealRoots(
    double a,
    double b,
    double c,
    double d,
    double tolerance)
{
    const double firstInvariant = b * b - 3 * a * c;
    const double secondInvariant = b * c - 9 * a * d;
    const double thirdInvariant = c * c - 3 * b * d;
    const double discriminant = secondInvariant * secondInvariant
        - 4 * firstInvariant * thirdInvariant;

    if (std::abs(firstInvariant) <= tolerance
        && std::abs(secondInvariant) <= tolerance) {
        const double root = -b / (3.0 * a);
        return { root, root, root };
    }
    if (std::abs(discriminant) <= tolerance) {
        const double ratio = secondInvariant / firstInvariant;
        const double firstRoot = -b / a + ratio;
        const double repeatedRoot = -ratio / 2.0;
        return { firstRoot, repeatedRoot, repeatedRoot };
    }
    if (discriminant < 0) {
        const double cosine = (2 * firstInvariant * b - 3 * a * secondInvariant)
            / (2 * firstInvariant * std::sqrt(firstInvariant));
        const double theta = std::acos(cosine);
        const double rootScale = std::sqrt(firstInvariant);
        return {
            (-b - 2 * rootScale * std::cos(theta / 3.0)) / (3 * a),
            (-b + rootScale * (std::cos(theta / 3.0)
                + std::sqrt(3.0) * std::sin(theta / 3.0))) / (3 * a),
            (-b + rootScale * (std::cos(theta / 3.0)
                - std::sqrt(3.0) * std::sin(theta / 3.0))) / (3 * a)
        };
    }

    const double firstTerm = firstInvariant * b
        + 1.5 * a * (-secondInvariant + std::sqrt(discriminant));
    const double secondTerm = firstInvariant * b
        + 1.5 * a * (-secondInvariant - std::sqrt(discriminant));
    return {
        (-b - std::cbrt(firstTerm) - std::cbrt(secondTerm)) / (3 * a)
    };
}

} // namespace ElasticityMath
