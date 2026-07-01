#include "math/MathUtils.h"

using namespace Eigen;

double vector_squareNorm(std::vector<Vector3d> vecs) {
	double norm = 0;
	for (const auto& v : vecs) {
		norm += v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	}
	return norm;
}

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

double f(double x, double a, double b, double c, double d) {
	double f = a * x * x * x + b * x * x + c * x + d;
	return f;
}

double df(double x, double a, double b, double c) {
	double df = 3 * a * x * x + 2 * b * x + c;
	return df;
}

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
		double v1 = f(pos[0], a, b, c, d);
		double v2 = f(pos[1], a, b, c, d);
		if ((v1) >= 0) {
			v1 = 0;
			results[1] = pos[0];
			results[2] = pos[0];
			//tempIndex = 2;
		}
		if ((v2) <= 0) {
			v2 = 0;
			results[1] = pos[1];
			results[2] = pos[1];
			//tempIndex = 2;
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
	double result[3];

	for (int i = 0; i < solves; i++) {
		double x1 = 0;
		int itCount = 0;
		do
		{
			if (itCount)
				x0 = x1;

			x1 = x0 - ((f(x0, a, b, c, d)) / (df(x0, a, b, c)));
			itCount++;

		} while (abs(x1 - x0) > EPS * EPS && itCount < 100000);
		results[tempIndex++] = x1;
		start = start + DX;
		x0 = start;
	}

	return results;
}

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
		double v1 = f(pos[0], a, b, c, d);
		double v2 = f(pos[1], a, b, c, d);
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
			if ((a < 0 && f(pos[0], a, b, c, d) > 0) || (a > 0 && f(pos[0], a, b, c, d) < 0)) {
				DX = -DX;
			}
		}
	}
	else if (delta == 0) {
		if (abs(f(specialPoint, a, b, c, d)) < EPS * EPS) {
			for (int i = 0; i < 3; i++) {
				double tempReuslt = specialPoint;
				results.push_back(tempReuslt);
			}
			return results;
		}
		if (a > 0) {
			if (f(specialPoint, a, b, c, d) > 0) {
				DX = 1;
			}
			else if (f(specialPoint, a, b, c, d) < 0) {
				DX = -1;
			}
		}
		else if (a < 0) {
			if (f(specialPoint, a, b, c, d) > 0) {
				DX = -1;
			}
			else if (f(specialPoint, a, b, c, d) < 0) {
				DX = 1;
			}
		}

	}

	double start = specialPoint - DX;
	double x0 = start;
	double result[3];

	for (int i = 0; i < solves; i++) {
		double x1 = 0;
		int itCount = 0;
		do
		{
			if (itCount)
				x0 = x1;

			x1 = x0 - ((f(x0, a, b, c, d)) / (df(x0, a, b, c)));
			itCount++;

		} while (abs(x1 - x0) > EPS * EPS && itCount < 10000);
		results.push_back(x1);
		start = start + DX;
		x0 = start;
	}

	return results;
}

std::vector<double> __SolverForCubicEquation(const double& a, const double& b, const double& c, const double& d, double EPS) {
	std::vector<double> resultsV;
	double A = b * b - 3 * a * c;
	double B = b * c - 9 * a * d;
	double C = c * c - 3 * b * d;
	double delta = B * B - 4 * A * C;
	double results[3];
	//num_solutions = 0;
	if (abs(A) == 0 && abs(B) == 0) {
		results[0] = -b / 3.0 / a;
		results[1] = results[0];
		results[2] = results[0];
		resultsV.push_back(results[0]);
		resultsV.push_back(results[1]);
		resultsV.push_back(results[2]);
		//num_solutions = 3;
	}
	else if (abs(delta) == 0) {
		double K = B / A;
		results[0] = -b / a + K;
		results[1] = -K / 2.0;
		results[2] = results[1];
		//num_solutions = 3;
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
		//num_solutions = 3;
		resultsV.push_back(results[0]);
		resultsV.push_back(results[1]);
		resultsV.push_back(results[2]);
	}
	else if (delta > 0) {
		double Y1 = A * b + 3 * a * (-B + sqrt(delta)) / 2;
		double Y2 = A * b + 3 * a * (-B - sqrt(delta)) / 2;
		//printf("Y1 = %f,  Y2 =  %f\n", Y1, Y2);
		results[0] = -b - cbrt(Y1) - cbrt(Y2);
		//num_solutions = 1;
		resultsV.push_back(results[0]);

	}
	return resultsV;
}

int linearIntersectTriangle(const Vector3d& mPoint0, const Vector3d& mPoint1, const Vector3d& mPoint2, const Vector3d& base, const Vector3d& direction, Vector3d& result)
{
	double u, v, tmp;
	Vector3d e1, e2, p, s, q;
	e1 = mPoint1 - mPoint0;
	e2 = mPoint2 - mPoint0;
	p = direction.cross(e2);
	tmp = p.dot(e1);
	//If the line is perpendicular to the normal of triangle.
	if (abs(tmp) < 1e-15)
	{
		p = e1.cross(e2);
		if (abs(p.dot(base - mPoint0)) < 1e-15)
			return 0;
		//The line is parallel to the plane.
		return 1;
	}
	s = base - mPoint0;
	u = p.dot(s) / tmp;
	if ((u < 0) || (u > 1))
		return 2;
	q = s.cross(e1);
	v = q.dot(direction) / tmp;
	if ((v < 0) || (v > 1) || (u + v > 1))
		return 2;

	result[0] = u;
	result[1] = v;
	result[2] = e2.dot(q) / tmp;
	return 3;
}

bool segmentHitTest(const Vector3d& point0, const Vector3d& point1, const Vector3d& point2, const Vector3d& seg_point0, const Vector3d& seg_point1)
{
	Vector3d ret;
	Vector3d direction = seg_point1 - seg_point0;
	int retFlag = linearIntersectTriangle(point0, point1, point2, seg_point0, direction, ret);
	if (retFlag != 3)
		return false;
	//check t. 0<= t <=1
	if ((ret[2] < 0) || (ret[2] > 1))
		return false;
	return true;
}