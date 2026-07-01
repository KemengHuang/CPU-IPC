#include "math/MathUtils.h"

using namespace Eigen;

// Sum of squared norms over a collection of 3-vectors.
double vector_squareNorm(const std::vector<Vector3d>& vecs) {
    double norm = 0;
    for (const auto& v : vecs) {
        norm += v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    }
    return norm;
}

// Dynamic-size column stacking.  Kept for non-hot callers.
MatrixXd vec_double(const MatrixXd& F) {
    const int cols = F.cols();
    const int rows = F.rows();
    const int nums = cols * rows;
    MatrixXd result(nums, 1);
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            result(i * rows + j, 0) = F(j, i);
        }
    }
    return result;
}

MatrixXf vec_float(const MatrixXf& F) {
    const int cols = F.cols();
    const int rows = F.rows();
    const int nums = cols * rows;
    MatrixXf result(nums, 1);
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            result(i * rows + j, 0) = F(j, i);
        }
    }
    return result;
}

// Fixed-size column stacking for the hot-path shapes.
Matrix<double, 9, 1> vec9(const Matrix3d& F) {
    Matrix<double, 9, 1> result;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result(i * 3 + j) = F(j, i);
        }
    }
    return result;
}

Matrix<double, 6, 1> vec6(const Matrix<double, 3, 2>& F) {
    Matrix<double, 6, 1> result;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            result(i * 3 + j) = F(j, i);
        }
    }
    return result;
}

Matrix<double, 12, 1> vec12(const Matrix<double, 3, 4>& F) {
    Matrix<double, 12, 1> result;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            result(i * 3 + j) = F(j, i);
        }
    }
    return result;
}

namespace {
    // Evaluate cubic a*x^3 + b*x^2 + c*x + d at x.
    double cubicValue(double x, double a, double b, double c, double d) {
        return a * x * x * x + b * x * x + c * x + d;
    }

    // Evaluate derivative 3*a*x^2 + 2*b*x + c at x.
    double cubicDerivative(double x, double a, double b, double c) {
        return 3 * a * x * x + 2 * b * x + c;
    }
}

// Newton-like solver used by the Stable Neo-Hookean eigen-projection.
// It estimates the number of real roots from the cubic discriminant and
// performs one or three Newton iterations starting from the inflection
// point and its surrounding extrema.
std::vector<double> NewtonSolverForCubicEquation_snk(const double& a, const double& b, const double& c, const double& d)
{
    double EPS = 1e-6;
    double DX = 0;
    std::vector<double> results(3);
    double specialPoint = -b / a / 3;
    double pos[2];
    int solves = 1;
    double delta = 4 * b * b - 12 * a * c;
    int tempIndex = 0;
    if (delta > 0) {
        pos[0] = (sqrt(delta) - 2 * b) / 6 / a;
        pos[1] = (-sqrt(delta) - 2 * b) / 6 / a;
        double v1 = cubicValue(pos[0], a, b, c, d);
        double v2 = cubicValue(pos[1], a, b, c, d);
        if ((v1) >= 0) {
            v1 = 0;
            results[1] = pos[0];
            results[2] = pos[0];
        }
        if ((v2) <= 0) {
            v2 = 0;
            results[1] = pos[1];
            results[2] = pos[1];
        }

        double sign = v1 * v2;
        DX = (pos[0] - pos[1]);
        if (sign < 0) {
            solves = 3;
        }
        else {
            if ((v2 <= 0)) {
                DX = -DX;
            }
        }
    }
    else {
        results[0] = specialPoint;
        results[1] = specialPoint;
        results[2] = specialPoint;
        return results;
    }

    double start = specialPoint - DX;
    double x0 = start;

    for (int i = 0; i < solves; i++) {
        double x1 = 0;
        int itCount = 0;
        do
        {
            if (itCount)
                x0 = x1;

            x1 = x0 - (cubicValue(x0, a, b, c, d) / cubicDerivative(x0, a, b, c));
            itCount++;

        } while (abs(x1 - x0) > EPS * EPS && itCount < 100000);
        results[tempIndex++] = x1;
        start = start + DX;
        x0 = start;
    }

    return results;
}

// General Newton-like cubic solver with explicit handling of the double-
// root case.  Returns all roots found; the number of roots may be one or
// three depending on the discriminant.
std::vector<double> NewtonSolverForCubicEquation(const double& a, const double& b, const double& c, const double& d)
{
    double EPS = 1e-6;
    double DX = 0;
    std::vector<double> results;
    double specialPoint = -b / a / 3;
    double pos[2];
    int solves = 1;
    double delta = 4 * b * b - 12 * a * c;
    if (delta > 0) {
        pos[0] = (sqrt(delta) - 2 * b) / 6 / a;
        pos[1] = (-sqrt(delta) - 2 * b) / 6 / a;
        double v1 = cubicValue(pos[0], a, b, c, d);
        double v2 = cubicValue(pos[1], a, b, c, d);
        if (abs(v1) < EPS * EPS) {
            v1 = 0;
        }
        if (abs(v2) < EPS * EPS) {
            v2 = 0;
        }
        double sign = v1 * v2;
        DX = (pos[0] - pos[1]);
        if (sign <= 0) {
            solves = 3;
        }
        else if (sign > 0) {
            if ((a < 0 && cubicValue(pos[0], a, b, c, d) > 0) || (a > 0 && cubicValue(pos[0], a, b, c, d) < 0)) {
                DX = -DX;
            }
        }
    }
    else if (delta == 0) {
        if (abs(cubicValue(specialPoint, a, b, c, d)) < EPS * EPS) {
            for (int i = 0; i < 3; i++) {
                double tempReuslt = specialPoint;
                results.push_back(tempReuslt);
            }
            return results;
        }
        if (a > 0) {
            if (cubicValue(specialPoint, a, b, c, d) > 0) {
                DX = 1;
            }
            else if (cubicValue(specialPoint, a, b, c, d) < 0) {
                DX = -1;
            }
        }
        else if (a < 0) {
            if (cubicValue(specialPoint, a, b, c, d) > 0) {
                DX = -1;
            }
            else if (cubicValue(specialPoint, a, b, c, d) < 0) {
                DX = 1;
            }
        }

    }

    double start = specialPoint - DX;
    double x0 = start;

    for (int i = 0; i < solves; i++) {
        double x1 = 0;
        int itCount = 0;
        do
        {
            if (itCount)
                x0 = x1;

            x1 = x0 - (cubicValue(x0, a, b, c, d) / cubicDerivative(x0, a, b, c));
            itCount++;

        } while (abs(x1 - x0) > EPS * EPS && itCount < 10000);
        results.push_back(x1);
        start = start + DX;
        x0 = start;
    }

    return results;
}

// Closed-form cubic solver based on the discriminant of the depressed
// cubic.  Returns one real root when delta > 0 and three (possibly
// repeated) real roots otherwise.
std::vector<double> __SolverForCubicEquation(const double& a, const double& b, const double& c, const double& d, double EPS) {
    std::vector<double> resultsV;
    double A = b * b - 3 * a * c;
    double B = b * c - 9 * a * d;
    double C = c * c - 3 * b * d;
    double delta = B * B - 4 * A * C;
    double results[3];
    if (abs(A) == 0 && abs(B) == 0) {
        results[0] = -b / 3.0 / a;
        results[1] = results[0];
        results[2] = results[0];
        resultsV.push_back(results[0]);
        resultsV.push_back(results[1]);
        resultsV.push_back(results[2]);
    }
    else if (abs(delta) == 0) {
        double K = B / A;
        results[0] = -b / a + K;
        results[1] = -K / 2.0;
        results[2] = results[1];
        resultsV.push_back(results[0]);
        resultsV.push_back(results[1]);
        resultsV.push_back(results[2]);
    }
    else if (delta < 0) {
        double T = (2 * A * b - 3 * a * B) / (2 * A * sqrt(A));
        double theta = acos(T);
        results[0] = (-b - 2 * sqrt(A) * cos(theta / 3.0)) / (3 * a);
        results[1] = (-b + sqrt(A) * (cos(theta / 3.0) + sqrt(3.0) * sin(theta / 3.0))) / (3 * a);
        results[2] = (-b + sqrt(A) * (cos(theta / 3.0) - sqrt(3.0) * sin(theta / 3.0))) / (3 * a);
        resultsV.push_back(results[0]);
        resultsV.push_back(results[1]);
        resultsV.push_back(results[2]);
    }
    else if (delta > 0) {
        double Y1 = A * b + 3 * a * (-B + sqrt(delta)) / 2;
        double Y2 = A * b + 3 * a * (-B - sqrt(delta)) / 2;
        results[0] = -b - cbrt(Y1) - cbrt(Y2);
        resultsV.push_back(results[0]);

    }
    return resultsV;
}

