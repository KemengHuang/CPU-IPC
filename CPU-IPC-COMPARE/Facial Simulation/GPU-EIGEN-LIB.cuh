#pragma once
#ifndef _GPU_EIGEN_LIB_CUH_
#define _GPU_EIGEN_LIB_CUH_
#include <cuda_runtime.h>
#include "math.h"
#include "Eigen_Data.h"

#define gone					1065353216
#define gsine_pi_over_eight		1053028117
#define gcosine_pi_over_eight   1064076127
#define gone_half				0.5f
#define gsmall_number			1.e-12f
#define gtiny_number			1.e-20f
#define gfour_gamma_squared		5.8284273147583007813f

union un { float f; unsigned int ui; };
namespace __GEIGEN__ {
	__device__ void __init_Mat3x3(Matrix3x3d& M, const double& val) {
		for (int i = 0;i < 3;i++) {
			for (int j = 0;j < 3;j++) {
				M.m[i][j] = val;
			}
		}
	}

	__device__ void __init_Mat9x9(Matrix9x9d& M, const double& val) {
		for (int i = 0;i < 9;i++) {
			for (int j = 0;j < 9;j++) {
				M.m[i][j] = val;
			}
		}
	}

	__device__ double __mabs(const double& a) {
		return a > 0 ? a : -a;
	}

	__device__ double __norm(const double3& n) {
		return sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
	}

	__device__ double3 __s_vec_multiply(const double3& a, double b) {
		return make_double3(a.x * b, a.y * b, a.z * b);
	}

	__device__ double3 __normalized(double3 n) {
		double norm = __norm(n);
		norm = 1 / norm;
		return __s_vec_multiply(n, norm);
	}

	__device__ double3 __add(double3 a, double3 b) {
		return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
	}

	__device__ Vector9 __add9(const Vector9& a, const Vector9& b) {
		Vector9 V;
		for (int i = 0;i < 9;i++) {
			V.v[i] = a.v[i] + b.v[i];
		}
		return V;
	}

	__device__ double3 __minus(double3 a, double3 b) {
		return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
	}

	__device__ double3 __v_vec_multiply(double3 a, double3 b) {
		return make_double3(a.x * b.x, a.y * b.y, a.z * b.z);
	}



	__device__ void __M_Mat_multiply(const Matrix3x3d& A, const Matrix3x3d& B, Matrix3x3d& output) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				double temp = 0;
				for (int k = 0; k < 3; k++) {
					temp += A.m[i][k] * B.m[k][j];
				}
				output.m[i][j] = temp;
			}
		}
	}

	__device__ double __Mat_Trace(const Matrix3x3d& A) {
		return A.m[0][0] + A.m[1][1] + A.m[2][2];
	}

	__device__ double3 __v_M_multiply(const double3& n, const Matrix3x3d& A) {
		double x = A.m[0][0] * n.x + A.m[1][0] * n.y + A.m[2][0] * n.z;
		double y = A.m[0][1] * n.x + A.m[1][1] * n.y + A.m[2][1] * n.z;
		double z = A.m[0][2] * n.x + A.m[1][2] * n.y + A.m[2][2] * n.z;
		return make_double3(x, y, z);
	}

	__device__ double3 __M_v_multiply(const Matrix3x3d& A, const double3& n) {
		double x = A.m[0][0] * n.x + A.m[0][1] * n.y + A.m[0][2] * n.z;
		double y = A.m[1][0] * n.x + A.m[1][1] * n.y + A.m[1][2] * n.z;
		double z = A.m[2][0] * n.x + A.m[2][1] * n.y + A.m[2][2] * n.z;
		return make_double3(x, y, z);
	}

	__device__ Vector12 __M12x9_v9_multiply(const Matrix12x9d& A, const Vector9& n) {
		Vector12 v12;
		for (int i = 0;i < 12;i++) {
			double temp = 0;
			for (int j = 0;j < 9;j++) {
				temp += A.m[i][j] * n.v[j];
			}
			v12.v[i] = temp;
		}
		return v12;
	}

	__device__ Vector12 __M12x12_v12_multiply(const Matrix12x12d& A, const Vector12& n) {
		Vector12 v12;
		for (int i = 0;i < 12;i++) {
			double temp = 0;
			for (int j = 0;j < 12;j++) {
				temp += A.m[i][j] * n.v[j];
			}
			v12.v[i] = temp;
		}
		return v12;
	}

	__device__ Matrix9x9d __S_Mat9x9_multiply(const Matrix9x9d& A, const double& B) {
		Matrix9x9d output;
		for (int i = 0;i < 9;i++) {
			for (int j = 0;j < 9;j++) {
				output.m[i][j] = A.m[i][j] * B;
			}
		}
		return output;
	}

	__device__ double __v_vec_dot(double3 a, double3 b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	__device__ Matrix3x3d __v_vec_toMat(double3 a, double3 b) {
		Matrix3x3d M;
		M.m[0][0] = a.x * b.x;M.m[0][1] = a.x * b.y;M.m[0][2] = a.x * b.z;
		M.m[1][0] = a.y * b.x;M.m[1][1] = a.y * b.y;M.m[1][2] = a.y * b.z;
		M.m[2][0] = a.z * b.x;M.m[2][1] = a.z * b.y;M.m[2][2] = a.z * b.z;
		return M;
	}

	__device__ Matrix9x9d __v9_vec9_toMat9x9(Vector9 a, Vector9 b) {
		Matrix9x9d M;
		for (int i = 0;i < 9;i++) {
			for (int j = 0;j < 9;j++) {
				M.m[i][j] = a.v[i] * b.v[j];
			}
		}
		return M;
	}

	__device__ Vector9 __s_vec9_multiply(Vector9 a, double b) {
		Vector9 V;
		for (int i = 0;i < 9;i++)
			V.v[i] = a.v[i] * b;
		return V;
	}

	__device__ void __Mat_add(const Matrix3x3d& A, const Matrix3x3d& B, Matrix3x3d& output) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				output.m[i][j] = A.m[i][j] + B.m[i][j];
	}

	__device__ Matrix9x9d __Mat9x9_add(const Matrix9x9d& A, const Matrix9x9d& B) {
		Matrix9x9d output;
		for (int i = 0; i < 9; i++)
			for (int j = 0; j < 9; j++)
				output.m[i][j] = A.m[i][j] + B.m[i][j];
		return output;
	}

	__device__ void __set_Mat_val(Matrix3x3d& M, const double& a00, const double& a01, const double& a02,
		const double& a10, const double& a11, const double& a12,
		const double& a20, const double& a21, const double& a22) {
		M.m[0][0] = a00;M.m[0][1] = a01;M.m[0][2] = a02;
		M.m[1][0] = a10;M.m[1][1] = a11;M.m[1][2] = a12;
		M.m[2][0] = a20;M.m[2][1] = a21;M.m[2][2] = a22;
	}

	__device__ void __init_Mat9x12_val(Matrix9x12d& M, const double& val) {
		for (int i = 0;i < 9;i++) {
			for (int j = 0;j < 12;j++) {
				M.m[i][j] = val;
			}
		}
	}

	__device__ Matrix3x3d __S_Mat_multiply(const Matrix3x3d& A, const double& B) {
		Matrix3x3d output;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				output.m[i][j] = A.m[i][j] * B;
		return output;
	}



	__device__ Matrix3x3d __Transpose3x3(Matrix3x3d input) {
		Matrix3x3d output;
		for (int i = 0;i < 3;i++) {
			for (int j = 0;j < 3;j++) {
				output.m[i][j] = input.m[j][i];
			}
		}
		return output;
		//output.m[0][0] = input.m[0][0];output.m[0][1] = input.m[1][0];output.m[0][2] = input.m[2][0];
		//output.m[1][0] = input.m[0][1];output.m[1][1] = input.m[1][1];output.m[1][2] = input.m[2][1];
		//output.m[2][0] = input.m[0][2];output.m[2][1] = input.m[1][2];output.m[2][2] = input.m[2][2];
	}

	__device__ Matrix12x9d __Transpose9x12(const Matrix9x12d& input) {
		Matrix12x9d output;
		for (int i = 0;i < 9;i++) {
			for (int j = 0;j < 12;j++) {
				output.m[j][i] = input.m[i][j];
			}
		}
		return output;
	}

	__device__ Matrix12x9d __M12x9_M9x9_Multiply(const Matrix12x9d& A, const Matrix9x9d& B) {
		Matrix12x9d output;
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < 9; j++) {
				double temp = 0;
				for (int k = 0; k < 9; k++) {
					temp += A.m[i][k] * B.m[k][j];
					
				}
				output.m[i][j] = temp;
			}
		}
		return output;
	}

	__device__ Matrix12x12d __M12x9_M9x12_Multiply(const Matrix12x9d& A, const Matrix9x12d& B) {
		Matrix12x12d output;
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < 12; j++) {
				double temp = 0;
				for (int k = 0; k < 9; k++) {
					temp += A.m[i][k] * B.m[k][j];
				}
				output.m[i][j] = temp;
			}
		}
		return output;
	}

	__device__ Matrix12x12d __s_M12x12_Multiply(const Matrix12x12d& A, const double& B) {
		Matrix12x12d output;
		for (int i = 0; i < 12; i++)
			for (int j = 0; j < 12; j++)
					output.m[i][j] = A.m[i][j] * B;
		return output;
	}

	__device__ void __Determiant(const Matrix3x3d& input, double& determinant) {
		determinant = input.m[0][0] * input.m[1][1] * input.m[2][2] +
			input.m[1][0] * input.m[2][1] * input.m[0][2] +
			input.m[2][0] * input.m[0][1] * input.m[1][2] -
			input.m[2][0] * input.m[1][1] * input.m[0][2] -
			input.m[0][0] * input.m[1][2] * input.m[2][1] -
			input.m[0][1] * input.m[1][0] * input.m[2][2];
	}

	__device__ void __Inverse(const Matrix3x3d& input, Matrix3x3d& output) {
		for (int i = 0;i < 3;i++) {
			for (int j = 0;j < 3;j++) {
				output.m[i][j] = input.m[i][j];
			}
		}

		int swapR[3], swapC[3];
		int pivot[3] = { 0 };

		for (int i = 0; i < 3; ++i) {
			int pr, pc;
			float maxValue = 0;
			for (int j = 0; j < 3; ++j) {
				if (pivot[j] != 1) {
					for (int k = 0; k < 3; ++k) {
						if (pivot[k] == 0) {
							if (output.m[j][k] > maxValue) {
								maxValue = output.m[j][k];
								pr = j;
								pc = k;
							}
						}
					}
				}
			}
			if (pr != pc) {
				float pv;
				for (int j = 0; j < 3; ++j) {
					pv = output.m[pr][j];
					output.m[pr][j] = output.m[pc][j];
					output.m[pc][j] = pv;
				}
			}
			swapC[i] = pc;
			swapR[i] = pr;
			++pivot[i];
			float inv = 1.f / output.m[pc][pc];
			output.m[pc][pc] = 1;
			for (int j = 0; j < 3; ++j) {
				output.m[pc][j] *= inv;
			}
			for (int j = 0; j < 3; ++j) {
				if (j != pc) {
					float powerRatio = output.m[j][pc];
					output.m[j][pc] = 0.f;
					for (int k = 0; k < 3; ++k) {
						output.m[j][k] -= output.m[pc][k] * powerRatio;
					}
				}
			}
		}
		for (int i = 0; i < 3; ++i) {
			if (swapR[i] != swapC[i]) {
				float pv;
				for (int j = 0; j < 3; ++j) {
					pv = output.m[j][swapC[i]];
					output.m[j][swapC[i]] = output.m[j][swapR[i]];
					output.m[j][swapR[i]] = pv;
				}
			}
		}
	}

	__device__ double __f(const double& x, const double& a, const double& b, const double& c, const double& d) {
		double f = a * x * x * x + b * x * x + c * x + d;
		return f;
	}

	__device__ double __df(const double& x, const double& a, const double& b, const double& c) {
		double df = 3 * a * x * x + 2 * b * x + c;
		return df;
	}

	__device__ void __NewtonSolverForCubicEquation(const double& a, const double& b, const double& c, const double& d, double* results, int& num_solutions)
	{
		double EPS = 1e-6;
		double DX = 0;
		//double results[3];
		num_solutions = 0;
		double specialPoint = -b / a / 3;
		double pos[2];
		int solves = 1;
		double delta = 4 * b * b - 12 * a * c;
		if (delta > 0) {
			pos[0] = (sqrt(delta) - 2 * b) / 6 / a;
			pos[1] = (-sqrt(delta) - 2 * b) / 6 / a;
			double v1 = __f(pos[0], a, b, c, d);
			double v2 = __f(pos[1], a, b, c, d);
			if (__mabs(v1) < 1e-10) {
				v1 = 0;
			}
			if (__mabs(v2) < 1e-10) {
				v2 = 0;
			}
			double sign = v1 * v2;
			DX = (pos[0] - pos[1]);
			if (sign <= 0) {
				solves = 3;
			}
			else if (sign > 0) {
				if ((a < 0 && __f(pos[0], a, b, c, d) > 0) || (a > 0 && __f(pos[0], a, b, c, d) < 0)) {
					DX = -DX;
				}
			}
		}
		else if (delta == 0) {
			if (__mabs(__f(specialPoint, a, b, c, d)) < 1e-10) {
				for (int i = 0; i < 3; i++) {
					double tempReuslt = specialPoint;
					results[num_solutions] = tempReuslt;
					num_solutions++;
				}
				return;
			}
			if (a > 0) {
				if (__f(specialPoint, a, b, c, d) > 0) {
					DX = 1;
				}
				else if (__f(specialPoint, a, b, c, d) < 0) {
					DX = -1;
				}
			}
			else if (a < 0) {
				if (__f(specialPoint, a, b, c, d) > 0) {
					DX = -1;
				}
				else if (__f(specialPoint, a, b, c, d) < 0) {
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

				x1 = x0 - ((__f(x0, a, b, c, d)) / (__df(x0, a, b, c)));
				itCount++;

			} while (__mabs(x1 - x0) > EPS && itCount < 1000);
			results[num_solutions] = (x1);
			num_solutions++;
			start = start + DX;
			x0 = start;
		}
	}

	__device__ __forceinline__
		void svd(
			float a11, float a12, float a13, float a21, float a22, float a23, float a31, float a32, float a33,			// input A     
			float& u11, float& u12, float& u13, float& u21, float& u22, float& u23, float& u31, float& u32, float& u33,	// output U      
			float& s11,
			//float &s12, float &s13, float &s21, 
			float& s22,
			//float &s23, float &s31, float &s32, 
			float& s33,	// output S
			float& v11, float& v12, float& v13, float& v21, float& v22, float& v23, float& v31, float& v32, float& v33	// output V
		)
	{
		un Sa11, Sa21, Sa31, Sa12, Sa22, Sa32, Sa13, Sa23, Sa33;
		un Su11, Su21, Su31, Su12, Su22, Su32, Su13, Su23, Su33;
		un Sv11, Sv21, Sv31, Sv12, Sv22, Sv32, Sv13, Sv23, Sv33;
		un Sc, Ss, Sch, Ssh;
		un Stmp1, Stmp2, Stmp3, Stmp4, Stmp5;
		un Ss11, Ss21, Ss31, Ss22, Ss32, Ss33;
		un Sqvs, Sqvvx, Sqvvy, Sqvvz;

		Sa11.f = a11;
		Sa12.f = a12;
		Sa13.f = a13;
		Sa21.f = a21;
		Sa22.f = a22;
		Sa23.f = a23;
		Sa31.f = a31;
		Sa32.f = a32;
		Sa33.f = a33;

		//###########################################################
		// Compute normal equations matrix
		//###########################################################

		Ss11.f = Sa11.f * Sa11.f;
		Stmp1.f = Sa21.f * Sa21.f;
		Ss11.f = __fadd_rn(Stmp1.f, Ss11.f);
		Stmp1.f = Sa31.f * Sa31.f;
		Ss11.f = __fadd_rn(Stmp1.f, Ss11.f);

		Ss21.f = Sa12.f * Sa11.f;
		Stmp1.f = Sa22.f * Sa21.f;
		Ss21.f = __fadd_rn(Stmp1.f, Ss21.f);
		Stmp1.f = Sa32.f * Sa31.f;
		Ss21.f = __fadd_rn(Stmp1.f, Ss21.f);

		Ss31.f = Sa13.f * Sa11.f;
		Stmp1.f = Sa23.f * Sa21.f;
		Ss31.f = __fadd_rn(Stmp1.f, Ss31.f);
		Stmp1.f = Sa33.f * Sa31.f;
		Ss31.f = __fadd_rn(Stmp1.f, Ss31.f);

		Ss22.f = Sa12.f * Sa12.f;
		Stmp1.f = Sa22.f * Sa22.f;
		Ss22.f = __fadd_rn(Stmp1.f, Ss22.f);
		Stmp1.f = Sa32.f * Sa32.f;
		Ss22.f = __fadd_rn(Stmp1.f, Ss22.f);

		Ss32.f = Sa13.f * Sa12.f;
		Stmp1.f = Sa23.f * Sa22.f;
		Ss32.f = __fadd_rn(Stmp1.f, Ss32.f);
		Stmp1.f = Sa33.f * Sa32.f;
		Ss32.f = __fadd_rn(Stmp1.f, Ss32.f);

		Ss33.f = Sa13.f * Sa13.f;
		Stmp1.f = Sa23.f * Sa23.f;
		Ss33.f = __fadd_rn(Stmp1.f, Ss33.f);
		Stmp1.f = Sa33.f * Sa33.f;
		Ss33.f = __fadd_rn(Stmp1.f, Ss33.f);

		Sqvs.f = 1.f;
		Sqvvx.f = 0.f;
		Sqvvy.f = 0.f;
		Sqvvz.f = 0.f;

		//###########################################################
		// Solve symmetric eigenproblem using Jacobi iteration
		//###########################################################
		for (int i = 0; i < 4; i++) {
			Ssh.f = Ss21.f * 0.5f;
			Stmp5.f = __fsub_rn(Ss11.f, Ss22.f);

			Stmp2.f = Ssh.f * Ssh.f;
			Stmp1.ui = (Stmp2.f >= gtiny_number) ? 0xffffffff : 0;
			Ssh.ui = Stmp1.ui & Ssh.ui;
			Sch.ui = Stmp1.ui & Stmp5.ui;
			Stmp2.ui = ~Stmp1.ui & gone;
			Sch.ui = Sch.ui | Stmp2.ui;

			Stmp1.f = Ssh.f * Ssh.f;
			Stmp2.f = Sch.f * Sch.f;
			Stmp3.f = __fadd_rn(Stmp1.f, Stmp2.f);
			Stmp4.f = __frsqrt_rn(Stmp3.f);

			Ssh.f = Stmp4.f * Ssh.f;
			Sch.f = Stmp4.f * Sch.f;
			Stmp1.f = gfour_gamma_squared * Stmp1.f;
			Stmp1.ui = (Stmp2.f <= Stmp1.f) ? 0xffffffff : 0;

			Stmp2.ui = gsine_pi_over_eight & Stmp1.ui;
			Ssh.ui = ~Stmp1.ui & Ssh.ui;
			Ssh.ui = Ssh.ui | Stmp2.ui;
			Stmp2.ui = gcosine_pi_over_eight & Stmp1.ui;
			Sch.ui = ~Stmp1.ui & Sch.ui;
			Sch.ui = Sch.ui | Stmp2.ui;

			Stmp1.f = Ssh.f * Ssh.f;
			Stmp2.f = Sch.f * Sch.f;
			Sc.f = __fsub_rn(Stmp2.f, Stmp1.f);
			Ss.f = Sch.f * Ssh.f;
			Ss.f = __fadd_rn(Ss.f, Ss.f);

#ifdef DEBUG_JACOBI_CONJUGATE
			printf("GPU s %.20g, c %.20g, sh %.20g, ch %.20g\n", Ss.f, Sc.f, Ssh.f, Sch.f);
#endif
			//###########################################################
			// Perform the actual Givens conjugation
			//###########################################################

			Stmp3.f = __fadd_rn(Stmp1.f, Stmp2.f);
			Ss33.f = Ss33.f * Stmp3.f;
			Ss31.f = Ss31.f * Stmp3.f;
			Ss32.f = Ss32.f * Stmp3.f;
			Ss33.f = Ss33.f * Stmp3.f;

			Stmp1.f = Ss.f * Ss31.f;
			Stmp2.f = Ss.f * Ss32.f;
			Ss31.f = Sc.f * Ss31.f;
			Ss32.f = Sc.f * Ss32.f;
			Ss31.f = __fadd_rn(Stmp2.f, Ss31.f);
			Ss32.f = __fsub_rn(Ss32.f, Stmp1.f);

			Stmp2.f = Ss.f * Ss.f;
			Stmp1.f = Ss22.f * Stmp2.f;
			Stmp3.f = Ss11.f * Stmp2.f;
			Stmp4.f = Sc.f * Sc.f;
			Ss11.f = Ss11.f * Stmp4.f;
			Ss22.f = Ss22.f * Stmp4.f;
			Ss11.f = __fadd_rn(Ss11.f, Stmp1.f);
			Ss22.f = __fadd_rn(Ss22.f, Stmp3.f);
			Stmp4.f = __fsub_rn(Stmp4.f, Stmp2.f);
			Stmp2.f = __fadd_rn(Ss21.f, Ss21.f);
			Ss21.f = Ss21.f * Stmp4.f;
			Stmp4.f = Sc.f * Ss.f;
			Stmp2.f = Stmp2.f * Stmp4.f;
			Stmp5.f = Stmp5.f * Stmp4.f;
			Ss11.f = __fadd_rn(Ss11.f, Stmp2.f);
			Ss21.f = __fsub_rn(Ss21.f, Stmp5.f);
			Ss22.f = __fsub_rn(Ss22.f, Stmp2.f);

#ifdef DEBUG_JACOBI_CONJUGATE
			printf("%.20g\n", Ss11.f);
			printf("%.20g %.20g\n", Ss21.f, Ss22.f);
			printf("%.20g %.20g %.20g\n", Ss31.f, Ss32.f, Ss33.f);
#endif

			//###########################################################
			// Compute the cumulative rotation, in quaternion form
			//###########################################################

			Stmp1.f = Ssh.f * Sqvvx.f;
			Stmp2.f = Ssh.f * Sqvvy.f;
			Stmp3.f = Ssh.f * Sqvvz.f;
			Ssh.f = Ssh.f * Sqvs.f;

			Sqvs.f = Sch.f * Sqvs.f;
			Sqvvx.f = Sch.f * Sqvvx.f;
			Sqvvy.f = Sch.f * Sqvvy.f;
			Sqvvz.f = Sch.f * Sqvvz.f;

			Sqvvz.f = __fadd_rn(Sqvvz.f, Ssh.f);
			Sqvs.f = __fsub_rn(Sqvs.f, Stmp3.f);
			Sqvvx.f = __fadd_rn(Sqvvx.f, Stmp2.f);
			Sqvvy.f = __fsub_rn(Sqvvy.f, Stmp1.f);

#ifdef DEBUG_JACOBI_CONJUGATE
			printf("GPU q %.20g %.20g %.20g %.20g\n", Sqvvx.f, Sqvvy.f, Sqvvz.f, Sqvs.f);
#endif

			//////////////////////////////////////////////////////////////////////////
			// (1->3)
			//////////////////////////////////////////////////////////////////////////
			Ssh.f = Ss32.f * 0.5f;
			Stmp5.f = __fsub_rn(Ss22.f, Ss33.f);

			Stmp2.f = Ssh.f * Ssh.f;
			Stmp1.ui = (Stmp2.f >= gtiny_number) ? 0xffffffff : 0;
			Ssh.ui = Stmp1.ui & Ssh.ui;
			Sch.ui = Stmp1.ui & Stmp5.ui;
			Stmp2.ui = ~Stmp1.ui & gone;
			Sch.ui = Sch.ui | Stmp2.ui;

			Stmp1.f = Ssh.f * Ssh.f;
			Stmp2.f = Sch.f * Sch.f;
			Stmp3.f = __fadd_rn(Stmp1.f, Stmp2.f);
			Stmp4.f = __frsqrt_rn(Stmp3.f);

			Ssh.f = Stmp4.f * Ssh.f;
			Sch.f = Stmp4.f * Sch.f;
			Stmp1.f = gfour_gamma_squared * Stmp1.f;
			Stmp1.ui = (Stmp2.f <= Stmp1.f) ? 0xffffffff : 0;

			Stmp2.ui = gsine_pi_over_eight & Stmp1.ui;
			Ssh.ui = ~Stmp1.ui & Ssh.ui;
			Ssh.ui = Ssh.ui | Stmp2.ui;
			Stmp2.ui = gcosine_pi_over_eight & Stmp1.ui;
			Sch.ui = ~Stmp1.ui & Sch.ui;
			Sch.ui = Sch.ui | Stmp2.ui;

			Stmp1.f = Ssh.f * Ssh.f;
			Stmp2.f = Sch.f * Sch.f;
			Sc.f = __fsub_rn(Stmp2.f, Stmp1.f);
			Ss.f = Sch.f * Ssh.f;
			Ss.f = __fadd_rn(Ss.f, Ss.f);

#ifdef DEBUG_JACOBI_CONJUGATE
			printf("GPU s %.20g, c %.20g, sh %.20g, ch %.20g\n", Ss.f, Sc.f, Ssh.f, Sch.f);
#endif

			//###########################################################
			// Perform the actual Givens conjugation
			//###########################################################

			Stmp3.f = __fadd_rn(Stmp1.f, Stmp2.f);
			Ss11.f = Ss11.f * Stmp3.f;
			Ss21.f = Ss21.f * Stmp3.f;
			Ss31.f = Ss31.f * Stmp3.f;
			Ss11.f = Ss11.f * Stmp3.f;

			Stmp1.f = Ss.f * Ss21.f;
			Stmp2.f = Ss.f * Ss31.f;
			Ss21.f = Sc.f * Ss21.f;
			Ss31.f = Sc.f * Ss31.f;
			Ss21.f = __fadd_rn(Stmp2.f, Ss21.f);
			Ss31.f = __fsub_rn(Ss31.f, Stmp1.f);

			Stmp2.f = Ss.f * Ss.f;
			Stmp1.f = Ss33.f * Stmp2.f;
			Stmp3.f = Ss22.f * Stmp2.f;
			Stmp4.f = Sc.f * Sc.f;
			Ss22.f = Ss22.f * Stmp4.f;
			Ss33.f = Ss33.f * Stmp4.f;
			Ss22.f = __fadd_rn(Ss22.f, Stmp1.f);
			Ss33.f = __fadd_rn(Ss33.f, Stmp3.f);
			Stmp4.f = __fsub_rn(Stmp4.f, Stmp2.f);
			Stmp2.f = __fadd_rn(Ss32.f, Ss32.f);
			Ss32.f = Ss32.f * Stmp4.f;
			Stmp4.f = Sc.f * Ss.f;
			Stmp2.f = Stmp2.f * Stmp4.f;
			Stmp5.f = Stmp5.f * Stmp4.f;
			Ss22.f = __fadd_rn(Ss22.f, Stmp2.f);
			Ss32.f = __fsub_rn(Ss32.f, Stmp5.f);
			Ss33.f = __fsub_rn(Ss33.f, Stmp2.f);

#ifdef DEBUG_JACOBI_CONJUGATE
			printf("%.20g\n", Ss11.f);
			printf("%.20g %.20g\n", Ss21.f, Ss22.f);
			printf("%.20g %.20g %.20g\n", Ss31.f, Ss32.f, Ss33.f);
#endif

			//###########################################################
			// Compute the cumulative rotation, in quaternion form
			//###########################################################

			Stmp1.f = Ssh.f * Sqvvx.f;
			Stmp2.f = Ssh.f * Sqvvy.f;
			Stmp3.f = Ssh.f * Sqvvz.f;
			Ssh.f = Ssh.f * Sqvs.f;

			Sqvs.f = Sch.f * Sqvs.f;
			Sqvvx.f = Sch.f * Sqvvx.f;
			Sqvvy.f = Sch.f * Sqvvy.f;
			Sqvvz.f = Sch.f * Sqvvz.f;

			Sqvvx.f = __fadd_rn(Sqvvx.f, Ssh.f);
			Sqvs.f = __fsub_rn(Sqvs.f, Stmp1.f);
			Sqvvy.f = __fadd_rn(Sqvvy.f, Stmp3.f);
			Sqvvz.f = __fsub_rn(Sqvvz.f, Stmp2.f);

#ifdef DEBUG_JACOBI_CONJUGATE
			printf("GPU q %.20g %.20g %.20g %.20g\n", Sqvvx.f, Sqvvy.f, Sqvvz.f, Sqvs.f);
#endif
#if 1
			//////////////////////////////////////////////////////////////////////////
			// 1 -> 2
			//////////////////////////////////////////////////////////////////////////

			Ssh.f = Ss31.f * 0.5f;
			Stmp5.f = __fsub_rn(Ss33.f, Ss11.f);

			Stmp2.f = Ssh.f * Ssh.f;
			Stmp1.ui = (Stmp2.f >= gtiny_number) ? 0xffffffff : 0;
			Ssh.ui = Stmp1.ui & Ssh.ui;
			Sch.ui = Stmp1.ui & Stmp5.ui;
			Stmp2.ui = ~Stmp1.ui & gone;
			Sch.ui = Sch.ui | Stmp2.ui;

			Stmp1.f = Ssh.f * Ssh.f;
			Stmp2.f = Sch.f * Sch.f;
			Stmp3.f = __fadd_rn(Stmp1.f, Stmp2.f);
			Stmp4.f = __frsqrt_rn(Stmp3.f);

			Ssh.f = Stmp4.f * Ssh.f;
			Sch.f = Stmp4.f * Sch.f;
			Stmp1.f = gfour_gamma_squared * Stmp1.f;
			Stmp1.ui = (Stmp2.f <= Stmp1.f) ? 0xffffffff : 0;

			Stmp2.ui = gsine_pi_over_eight & Stmp1.ui;
			Ssh.ui = ~Stmp1.ui & Ssh.ui;
			Ssh.ui = Ssh.ui | Stmp2.ui;
			Stmp2.ui = gcosine_pi_over_eight & Stmp1.ui;
			Sch.ui = ~Stmp1.ui & Sch.ui;
			Sch.ui = Sch.ui | Stmp2.ui;

			Stmp1.f = Ssh.f * Ssh.f;
			Stmp2.f = Sch.f * Sch.f;
			Sc.f = __fsub_rn(Stmp2.f, Stmp1.f);
			Ss.f = Sch.f * Ssh.f;
			Ss.f = __fadd_rn(Ss.f, Ss.f);

#ifdef DEBUG_JACOBI_CONJUGATE
			printf("GPU s %.20g, c %.20g, sh %.20g, ch %.20g\n", Ss.f, Sc.f, Ssh.f, Sch.f);
#endif

			//###########################################################
			// Perform the actual Givens conjugation
			//###########################################################

			Stmp3.f = __fadd_rn(Stmp1.f, Stmp2.f);
			Ss22.f = Ss22.f * Stmp3.f;
			Ss32.f = Ss32.f * Stmp3.f;
			Ss21.f = Ss21.f * Stmp3.f;
			Ss22.f = Ss22.f * Stmp3.f;

			Stmp1.f = Ss.f * Ss32.f;
			Stmp2.f = Ss.f * Ss21.f;
			Ss32.f = Sc.f * Ss32.f;
			Ss21.f = Sc.f * Ss21.f;
			Ss32.f = __fadd_rn(Stmp2.f, Ss32.f);
			Ss21.f = __fsub_rn(Ss21.f, Stmp1.f);

			Stmp2.f = Ss.f * Ss.f;
			Stmp1.f = Ss11.f * Stmp2.f;
			Stmp3.f = Ss33.f * Stmp2.f;
			Stmp4.f = Sc.f * Sc.f;
			Ss33.f = Ss33.f * Stmp4.f;
			Ss11.f = Ss11.f * Stmp4.f;
			Ss33.f = __fadd_rn(Ss33.f, Stmp1.f);
			Ss11.f = __fadd_rn(Ss11.f, Stmp3.f);
			Stmp4.f = __fsub_rn(Stmp4.f, Stmp2.f);
			Stmp2.f = __fadd_rn(Ss31.f, Ss31.f);
			Ss31.f = Ss31.f * Stmp4.f;
			Stmp4.f = Sc.f * Ss.f;
			Stmp2.f = Stmp2.f * Stmp4.f;
			Stmp5.f = Stmp5.f * Stmp4.f;
			Ss33.f = __fadd_rn(Ss33.f, Stmp2.f);
			Ss31.f = __fsub_rn(Ss31.f, Stmp5.f);
			Ss11.f = __fsub_rn(Ss11.f, Stmp2.f);

#ifdef DEBUG_JACOBI_CONJUGATE
			printf("%.20g\n", Ss11.f);
			printf("%.20g %.20g\n", Ss21.f, Ss22.f);
			printf("%.20g %.20g %.20g\n", Ss31.f, Ss32.f, Ss33.f);
#endif

			//###########################################################
			// Compute the cumulative rotation, in quaternion form
			//###########################################################

			Stmp1.f = Ssh.f * Sqvvx.f;
			Stmp2.f = Ssh.f * Sqvvy.f;
			Stmp3.f = Ssh.f * Sqvvz.f;
			Ssh.f = Ssh.f * Sqvs.f;

			Sqvs.f = Sch.f * Sqvs.f;
			Sqvvx.f = Sch.f * Sqvvx.f;
			Sqvvy.f = Sch.f * Sqvvy.f;
			Sqvvz.f = Sch.f * Sqvvz.f;

			Sqvvy.f = __fadd_rn(Sqvvy.f, Ssh.f);
			Sqvs.f = __fsub_rn(Sqvs.f, Stmp2.f);
			Sqvvz.f = __fadd_rn(Sqvvz.f, Stmp1.f);
			Sqvvx.f = __fsub_rn(Sqvvx.f, Stmp3.f);
#endif
		}

		//###########################################################
		// Normalize quaternion for matrix V
		//###########################################################

		Stmp2.f = Sqvs.f * Sqvs.f;
		Stmp1.f = Sqvvx.f * Sqvvx.f;
		Stmp2.f = __fadd_rn(Stmp1.f, Stmp2.f);
		Stmp1.f = Sqvvy.f * Sqvvy.f;
		Stmp2.f = __fadd_rn(Stmp1.f, Stmp2.f);
		Stmp1.f = Sqvvz.f * Sqvvz.f;
		Stmp2.f = __fadd_rn(Stmp1.f, Stmp2.f);

		Stmp1.f = __frsqrt_rn(Stmp2.f);
		Stmp4.f = Stmp1.f * 0.5f;
		Stmp3.f = Stmp1.f * Stmp4.f;
		Stmp3.f = Stmp1.f * Stmp3.f;
		Stmp3.f = Stmp2.f * Stmp3.f;
		Stmp1.f = __fadd_rn(Stmp1.f, Stmp4.f);
		Stmp1.f = __fsub_rn(Stmp1.f, Stmp3.f);

		Sqvs.f = Sqvs.f * Stmp1.f;
		Sqvvx.f = Sqvvx.f * Stmp1.f;
		Sqvvy.f = Sqvvy.f * Stmp1.f;
		Sqvvz.f = Sqvvz.f * Stmp1.f;

		//###########################################################
		// Transform quaternion to matrix V
		//###########################################################

		Stmp1.f = Sqvvx.f * Sqvvx.f;
		Stmp2.f = Sqvvy.f * Sqvvy.f;
		Stmp3.f = Sqvvz.f * Sqvvz.f;
		Sv11.f = Sqvs.f * Sqvs.f;
		Sv22.f = __fsub_rn(Sv11.f, Stmp1.f);
		Sv33.f = __fsub_rn(Sv22.f, Stmp2.f);
		Sv33.f = __fadd_rn(Sv33.f, Stmp3.f);
		Sv22.f = __fadd_rn(Sv22.f, Stmp2.f);
		Sv22.f = __fsub_rn(Sv22.f, Stmp3.f);
		Sv11.f = __fadd_rn(Sv11.f, Stmp1.f);
		Sv11.f = __fsub_rn(Sv11.f, Stmp2.f);
		Sv11.f = __fsub_rn(Sv11.f, Stmp3.f);
		Stmp1.f = __fadd_rn(Sqvvx.f, Sqvvx.f);
		Stmp2.f = __fadd_rn(Sqvvy.f, Sqvvy.f);
		Stmp3.f = __fadd_rn(Sqvvz.f, Sqvvz.f);
		Sv32.f = Sqvs.f * Stmp1.f;
		Sv13.f = Sqvs.f * Stmp2.f;
		Sv21.f = Sqvs.f * Stmp3.f;
		Stmp1.f = Sqvvy.f * Stmp1.f;
		Stmp2.f = Sqvvz.f * Stmp2.f;
		Stmp3.f = Sqvvx.f * Stmp3.f;
		Sv12.f = __fsub_rn(Stmp1.f, Sv21.f);
		Sv23.f = __fsub_rn(Stmp2.f, Sv32.f);
		Sv31.f = __fsub_rn(Stmp3.f, Sv13.f);
		Sv21.f = __fadd_rn(Stmp1.f, Sv21.f);
		Sv32.f = __fadd_rn(Stmp2.f, Sv32.f);
		Sv13.f = __fadd_rn(Stmp3.f, Sv13.f);

		///###########################################################
		// Multiply (from the right) with V
		//###########################################################

		Stmp2.f = Sa12.f;
		Stmp3.f = Sa13.f;
		Sa12.f = Sv12.f * Sa11.f;
		Sa13.f = Sv13.f * Sa11.f;
		Sa11.f = Sv11.f * Sa11.f;
		Stmp1.f = Sv21.f * Stmp2.f;
		Sa11.f = __fadd_rn(Sa11.f, Stmp1.f);
		Stmp1.f = Sv31.f * Stmp3.f;
		Sa11.f = __fadd_rn(Sa11.f, Stmp1.f);
		Stmp1.f = Sv22.f * Stmp2.f;
		Sa12.f = __fadd_rn(Sa12.f, Stmp1.f);
		Stmp1.f = Sv32.f * Stmp3.f;
		Sa12.f = __fadd_rn(Sa12.f, Stmp1.f);
		Stmp1.f = Sv23.f * Stmp2.f;
		Sa13.f = __fadd_rn(Sa13.f, Stmp1.f);
		Stmp1.f = Sv33.f * Stmp3.f;
		Sa13.f = __fadd_rn(Sa13.f, Stmp1.f);

		Stmp2.f = Sa22.f;
		Stmp3.f = Sa23.f;
		Sa22.f = Sv12.f * Sa21.f;
		Sa23.f = Sv13.f * Sa21.f;
		Sa21.f = Sv11.f * Sa21.f;
		Stmp1.f = Sv21.f * Stmp2.f;
		Sa21.f = __fadd_rn(Sa21.f, Stmp1.f);
		Stmp1.f = Sv31.f * Stmp3.f;
		Sa21.f = __fadd_rn(Sa21.f, Stmp1.f);
		Stmp1.f = Sv22.f * Stmp2.f;
		Sa22.f = __fadd_rn(Sa22.f, Stmp1.f);
		Stmp1.f = Sv32.f * Stmp3.f;
		Sa22.f = __fadd_rn(Sa22.f, Stmp1.f);
		Stmp1.f = Sv23.f * Stmp2.f;
		Sa23.f = __fadd_rn(Sa23.f, Stmp1.f);
		Stmp1.f = Sv33.f * Stmp3.f;
		Sa23.f = __fadd_rn(Sa23.f, Stmp1.f);

		Stmp2.f = Sa32.f;
		Stmp3.f = Sa33.f;
		Sa32.f = Sv12.f * Sa31.f;
		Sa33.f = Sv13.f * Sa31.f;
		Sa31.f = Sv11.f * Sa31.f;
		Stmp1.f = Sv21.f * Stmp2.f;
		Sa31.f = __fadd_rn(Sa31.f, Stmp1.f);
		Stmp1.f = Sv31.f * Stmp3.f;
		Sa31.f = __fadd_rn(Sa31.f, Stmp1.f);
		Stmp1.f = Sv22.f * Stmp2.f;
		Sa32.f = __fadd_rn(Sa32.f, Stmp1.f);
		Stmp1.f = Sv32.f * Stmp3.f;
		Sa32.f = __fadd_rn(Sa32.f, Stmp1.f);
		Stmp1.f = Sv23.f * Stmp2.f;
		Sa33.f = __fadd_rn(Sa33.f, Stmp1.f);
		Stmp1.f = Sv33.f * Stmp3.f;
		Sa33.f = __fadd_rn(Sa33.f, Stmp1.f);

		//###########################################################
		// Permute columns such that the singular values are sorted
		//###########################################################

		Stmp1.f = Sa11.f * Sa11.f;
		Stmp4.f = Sa21.f * Sa21.f;
		Stmp1.f = __fadd_rn(Stmp1.f, Stmp4.f);
		Stmp4.f = Sa31.f * Sa31.f;
		Stmp1.f = __fadd_rn(Stmp1.f, Stmp4.f);

		Stmp2.f = Sa12.f * Sa12.f;
		Stmp4.f = Sa22.f * Sa22.f;
		Stmp2.f = __fadd_rn(Stmp2.f, Stmp4.f);
		Stmp4.f = Sa32.f * Sa32.f;
		Stmp2.f = __fadd_rn(Stmp2.f, Stmp4.f);

		Stmp3.f = Sa13.f * Sa13.f;
		Stmp4.f = Sa23.f * Sa23.f;
		Stmp3.f = __fadd_rn(Stmp3.f, Stmp4.f);
		Stmp4.f = Sa33.f * Sa33.f;
		Stmp3.f = __fadd_rn(Stmp3.f, Stmp4.f);

		// Swap columns 1-2 if necessary

		Stmp4.ui = (Stmp1.f < Stmp2.f) ? 0xffffffff : 0;
		Stmp5.ui = Sa11.ui ^ Sa12.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sa11.ui = Sa11.ui ^ Stmp5.ui;
		Sa12.ui = Sa12.ui ^ Stmp5.ui;

		Stmp5.ui = Sa21.ui ^ Sa22.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sa21.ui = Sa21.ui ^ Stmp5.ui;
		Sa22.ui = Sa22.ui ^ Stmp5.ui;

		Stmp5.ui = Sa31.ui ^ Sa32.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sa31.ui = Sa31.ui ^ Stmp5.ui;
		Sa32.ui = Sa32.ui ^ Stmp5.ui;

		Stmp5.ui = Sv11.ui ^ Sv12.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sv11.ui = Sv11.ui ^ Stmp5.ui;
		Sv12.ui = Sv12.ui ^ Stmp5.ui;

		Stmp5.ui = Sv21.ui ^ Sv22.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sv21.ui = Sv21.ui ^ Stmp5.ui;
		Sv22.ui = Sv22.ui ^ Stmp5.ui;

		Stmp5.ui = Sv31.ui ^ Sv32.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sv31.ui = Sv31.ui ^ Stmp5.ui;
		Sv32.ui = Sv32.ui ^ Stmp5.ui;

		Stmp5.ui = Stmp1.ui ^ Stmp2.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Stmp1.ui = Stmp1.ui ^ Stmp5.ui;
		Stmp2.ui = Stmp2.ui ^ Stmp5.ui;

		// If columns 1-2 have been swapped, negate 2nd column of A and V so that V is still a rotation

		Stmp5.f = -2.f;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Stmp4.f = 1.f;
		Stmp4.f = __fadd_rn(Stmp4.f, Stmp5.f);

		Sa12.f = Sa12.f * Stmp4.f;
		Sa22.f = Sa22.f * Stmp4.f;
		Sa32.f = Sa32.f * Stmp4.f;

		Sv12.f = Sv12.f * Stmp4.f;
		Sv22.f = Sv22.f * Stmp4.f;
		Sv32.f = Sv32.f * Stmp4.f;

		// Swap columns 1-3 if necessary

		Stmp4.ui = (Stmp1.f < Stmp3.f) ? 0xffffffff : 0;
		Stmp5.ui = Sa11.ui ^ Sa13.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sa11.ui = Sa11.ui ^ Stmp5.ui;
		Sa13.ui = Sa13.ui ^ Stmp5.ui;

		Stmp5.ui = Sa21.ui ^ Sa23.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sa21.ui = Sa21.ui ^ Stmp5.ui;
		Sa23.ui = Sa23.ui ^ Stmp5.ui;

		Stmp5.ui = Sa31.ui ^ Sa33.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sa31.ui = Sa31.ui ^ Stmp5.ui;
		Sa33.ui = Sa33.ui ^ Stmp5.ui;

		Stmp5.ui = Sv11.ui ^ Sv13.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sv11.ui = Sv11.ui ^ Stmp5.ui;
		Sv13.ui = Sv13.ui ^ Stmp5.ui;

		Stmp5.ui = Sv21.ui ^ Sv23.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sv21.ui = Sv21.ui ^ Stmp5.ui;
		Sv23.ui = Sv23.ui ^ Stmp5.ui;

		Stmp5.ui = Sv31.ui ^ Sv33.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sv31.ui = Sv31.ui ^ Stmp5.ui;
		Sv33.ui = Sv33.ui ^ Stmp5.ui;

		Stmp5.ui = Stmp1.ui ^ Stmp3.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Stmp1.ui = Stmp1.ui ^ Stmp5.ui;
		Stmp3.ui = Stmp3.ui ^ Stmp5.ui;

		// If columns 1-3 have been swapped, negate 1st column of A and V so that V is still a rotation

		Stmp5.f = -2.f;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Stmp4.f = 1.f;
		Stmp4.f = __fadd_rn(Stmp4.f, Stmp5.f);

		Sa11.f = Sa11.f * Stmp4.f;
		Sa21.f = Sa21.f * Stmp4.f;
		Sa31.f = Sa31.f * Stmp4.f;

		Sv11.f = Sv11.f * Stmp4.f;
		Sv21.f = Sv21.f * Stmp4.f;
		Sv31.f = Sv31.f * Stmp4.f;

		// Swap columns 2-3 if necessary

		Stmp4.ui = (Stmp2.f < Stmp3.f) ? 0xffffffff : 0;
		Stmp5.ui = Sa12.ui ^ Sa13.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sa12.ui = Sa12.ui ^ Stmp5.ui;
		Sa13.ui = Sa13.ui ^ Stmp5.ui;

		Stmp5.ui = Sa22.ui ^ Sa23.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sa22.ui = Sa22.ui ^ Stmp5.ui;
		Sa23.ui = Sa23.ui ^ Stmp5.ui;

		Stmp5.ui = Sa32.ui ^ Sa33.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sa32.ui = Sa32.ui ^ Stmp5.ui;
		Sa33.ui = Sa33.ui ^ Stmp5.ui;

		Stmp5.ui = Sv12.ui ^ Sv13.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sv12.ui = Sv12.ui ^ Stmp5.ui;
		Sv13.ui = Sv13.ui ^ Stmp5.ui;

		Stmp5.ui = Sv22.ui ^ Sv23.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sv22.ui = Sv22.ui ^ Stmp5.ui;
		Sv23.ui = Sv23.ui ^ Stmp5.ui;

		Stmp5.ui = Sv32.ui ^ Sv33.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Sv32.ui = Sv32.ui ^ Stmp5.ui;
		Sv33.ui = Sv33.ui ^ Stmp5.ui;

		Stmp5.ui = Stmp2.ui ^ Stmp3.ui;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Stmp2.ui = Stmp2.ui ^ Stmp5.ui;
		Stmp3.ui = Stmp3.ui ^ Stmp5.ui;

		// If columns 2-3 have been swapped, negate 3rd column of A and V so that V is still a rotation

		Stmp5.f = -2.f;
		Stmp5.ui = Stmp5.ui & Stmp4.ui;
		Stmp4.f = 1.f;
		Stmp4.f = __fadd_rn(Stmp4.f, Stmp5.f);

		Sa13.f = Sa13.f * Stmp4.f;
		Sa23.f = Sa23.f * Stmp4.f;
		Sa33.f = Sa33.f * Stmp4.f;

		Sv13.f = Sv13.f * Stmp4.f;
		Sv23.f = Sv23.f * Stmp4.f;
		Sv33.f = Sv33.f * Stmp4.f;

		//###########################################################
		// Construct QR factorization of A*V (=U*D) using Givens rotations
		//###########################################################

		Su11.f = 1.f;
		Su12.f = 0.f;
		Su13.f = 0.f;
		Su21.f = 0.f;
		Su22.f = 1.f;
		Su23.f = 0.f;
		Su31.f = 0.f;
		Su32.f = 0.f;
		Su33.f = 1.f;

		Ssh.f = Sa21.f * Sa21.f;
		Ssh.ui = (Ssh.f >= gsmall_number) ? 0xffffffff : 0;
		Ssh.ui = Ssh.ui & Sa21.ui;

		Stmp5.f = 0.f;
		Sch.f = __fsub_rn(Stmp5.f, Sa11.f);
		Sch.f = max(Sch.f, Sa11.f);
		Sch.f = max(Sch.f, gsmall_number);
		Stmp5.ui = (Sa11.f >= Stmp5.f) ? 0xffffffff : 0;

		Stmp1.f = Sch.f * Sch.f;
		Stmp2.f = Ssh.f * Ssh.f;
		Stmp2.f = __fadd_rn(Stmp1.f, Stmp2.f);
		Stmp1.f = __frsqrt_rn(Stmp2.f);

		Stmp4.f = Stmp1.f * 0.5f;
		Stmp3.f = Stmp1.f * Stmp4.f;
		Stmp3.f = Stmp1.f * Stmp3.f;
		Stmp3.f = Stmp2.f * Stmp3.f;
		Stmp1.f = __fadd_rn(Stmp1.f, Stmp4.f);
		Stmp1.f = __fsub_rn(Stmp1.f, Stmp3.f);
		Stmp1.f = Stmp1.f * Stmp2.f;

		Sch.f = __fadd_rn(Sch.f, Stmp1.f);

		Stmp1.ui = ~Stmp5.ui & Ssh.ui;
		Stmp2.ui = ~Stmp5.ui & Sch.ui;
		Sch.ui = Stmp5.ui & Sch.ui;
		Ssh.ui = Stmp5.ui & Ssh.ui;
		Sch.ui = Sch.ui | Stmp1.ui;
		Ssh.ui = Ssh.ui | Stmp2.ui;

		Stmp1.f = Sch.f * Sch.f;
		Stmp2.f = Ssh.f * Ssh.f;
		Stmp2.f = __fadd_rn(Stmp1.f, Stmp2.f);
		Stmp1.f = __frsqrt_rn(Stmp2.f);

		Stmp4.f = Stmp1.f * 0.5f;
		Stmp3.f = Stmp1.f * Stmp4.f;
		Stmp3.f = Stmp1.f * Stmp3.f;
		Stmp3.f = Stmp2.f * Stmp3.f;
		Stmp1.f = __fadd_rn(Stmp1.f, Stmp4.f);
		Stmp1.f = __fsub_rn(Stmp1.f, Stmp3.f);

		Sch.f = Sch.f * Stmp1.f;
		Ssh.f = Ssh.f * Stmp1.f;

		Sc.f = Sch.f * Sch.f;
		Ss.f = Ssh.f * Ssh.f;
		Sc.f = __fsub_rn(Sc.f, Ss.f);
		Ss.f = Ssh.f * Sch.f;
		Ss.f = __fadd_rn(Ss.f, Ss.f);

		//###########################################################
		// Rotate matrix A
		//###########################################################

		Stmp1.f = Ss.f * Sa11.f;
		Stmp2.f = Ss.f * Sa21.f;
		Sa11.f = Sc.f * Sa11.f;
		Sa21.f = Sc.f * Sa21.f;
		Sa11.f = __fadd_rn(Sa11.f, Stmp2.f);
		Sa21.f = __fsub_rn(Sa21.f, Stmp1.f);

		Stmp1.f = Ss.f * Sa12.f;
		Stmp2.f = Ss.f * Sa22.f;
		Sa12.f = Sc.f * Sa12.f;
		Sa22.f = Sc.f * Sa22.f;
		Sa12.f = __fadd_rn(Sa12.f, Stmp2.f);
		Sa22.f = __fsub_rn(Sa22.f, Stmp1.f);

		Stmp1.f = Ss.f * Sa13.f;
		Stmp2.f = Ss.f * Sa23.f;
		Sa13.f = Sc.f * Sa13.f;
		Sa23.f = Sc.f * Sa23.f;
		Sa13.f = __fadd_rn(Sa13.f, Stmp2.f);
		Sa23.f = __fsub_rn(Sa23.f, Stmp1.f);

		//###########################################################
		// Update matrix U
		//###########################################################

		Stmp1.f = Ss.f * Su11.f;
		Stmp2.f = Ss.f * Su12.f;
		Su11.f = Sc.f * Su11.f;
		Su12.f = Sc.f * Su12.f;
		Su11.f = __fadd_rn(Su11.f, Stmp2.f);
		Su12.f = __fsub_rn(Su12.f, Stmp1.f);

		Stmp1.f = Ss.f * Su21.f;
		Stmp2.f = Ss.f * Su22.f;
		Su21.f = Sc.f * Su21.f;
		Su22.f = Sc.f * Su22.f;
		Su21.f = __fadd_rn(Su21.f, Stmp2.f);
		Su22.f = __fsub_rn(Su22.f, Stmp1.f);

		Stmp1.f = Ss.f * Su31.f;
		Stmp2.f = Ss.f * Su32.f;
		Su31.f = Sc.f * Su31.f;
		Su32.f = Sc.f * Su32.f;
		Su31.f = __fadd_rn(Su31.f, Stmp2.f);
		Su32.f = __fsub_rn(Su32.f, Stmp1.f);

		// Second Givens rotation

		Ssh.f = Sa31.f * Sa31.f;
		Ssh.ui = (Ssh.f >= gsmall_number) ? 0xffffffff : 0;
		Ssh.ui = Ssh.ui & Sa31.ui;

		Stmp5.f = 0.f;
		Sch.f = __fsub_rn(Stmp5.f, Sa11.f);
		Sch.f = max(Sch.f, Sa11.f);
		Sch.f = max(Sch.f, gsmall_number);
		Stmp5.ui = (Sa11.f >= Stmp5.f) ? 0xffffffff : 0;

		Stmp1.f = Sch.f * Sch.f;
		Stmp2.f = Ssh.f * Ssh.f;
		Stmp2.f = __fadd_rn(Stmp1.f, Stmp2.f);
		Stmp1.f = __frsqrt_rn(Stmp2.f);

		Stmp4.f = Stmp1.f * 0.5;
		Stmp3.f = Stmp1.f * Stmp4.f;
		Stmp3.f = Stmp1.f * Stmp3.f;
		Stmp3.f = Stmp2.f * Stmp3.f;
		Stmp1.f = __fadd_rn(Stmp1.f, Stmp4.f);
		Stmp1.f = __fsub_rn(Stmp1.f, Stmp3.f);
		Stmp1.f = Stmp1.f * Stmp2.f;

		Sch.f = __fadd_rn(Sch.f, Stmp1.f);

		Stmp1.ui = ~Stmp5.ui & Ssh.ui;
		Stmp2.ui = ~Stmp5.ui & Sch.ui;
		Sch.ui = Stmp5.ui & Sch.ui;
		Ssh.ui = Stmp5.ui & Ssh.ui;
		Sch.ui = Sch.ui | Stmp1.ui;
		Ssh.ui = Ssh.ui | Stmp2.ui;

		Stmp1.f = Sch.f * Sch.f;
		Stmp2.f = Ssh.f * Ssh.f;
		Stmp2.f = __fadd_rn(Stmp1.f, Stmp2.f);
		Stmp1.f = __frsqrt_rn(Stmp2.f);

		Stmp4.f = Stmp1.f * 0.5f;
		Stmp3.f = Stmp1.f * Stmp4.f;
		Stmp3.f = Stmp1.f * Stmp3.f;
		Stmp3.f = Stmp2.f * Stmp3.f;
		Stmp1.f = __fadd_rn(Stmp1.f, Stmp4.f);
		Stmp1.f = __fsub_rn(Stmp1.f, Stmp3.f);

		Sch.f = Sch.f * Stmp1.f;
		Ssh.f = Ssh.f * Stmp1.f;

		Sc.f = Sch.f * Sch.f;
		Ss.f = Ssh.f * Ssh.f;
		Sc.f = __fsub_rn(Sc.f, Ss.f);
		Ss.f = Ssh.f * Sch.f;
		Ss.f = __fadd_rn(Ss.f, Ss.f);

		//###########################################################
		// Rotate matrix A
		//###########################################################

		Stmp1.f = Ss.f * Sa11.f;
		Stmp2.f = Ss.f * Sa31.f;
		Sa11.f = Sc.f * Sa11.f;
		Sa31.f = Sc.f * Sa31.f;
		Sa11.f = __fadd_rn(Sa11.f, Stmp2.f);
		Sa31.f = __fsub_rn(Sa31.f, Stmp1.f);

		Stmp1.f = Ss.f * Sa12.f;
		Stmp2.f = Ss.f * Sa32.f;
		Sa12.f = Sc.f * Sa12.f;
		Sa32.f = Sc.f * Sa32.f;
		Sa12.f = __fadd_rn(Sa12.f, Stmp2.f);
		Sa32.f = __fsub_rn(Sa32.f, Stmp1.f);

		Stmp1.f = Ss.f * Sa13.f;
		Stmp2.f = Ss.f * Sa33.f;
		Sa13.f = Sc.f * Sa13.f;
		Sa33.f = Sc.f * Sa33.f;
		Sa13.f = __fadd_rn(Sa13.f, Stmp2.f);
		Sa33.f = __fsub_rn(Sa33.f, Stmp1.f);

		//###########################################################
		// Update matrix U
		//###########################################################

		Stmp1.f = Ss.f * Su11.f;
		Stmp2.f = Ss.f * Su13.f;
		Su11.f = Sc.f * Su11.f;
		Su13.f = Sc.f * Su13.f;
		Su11.f = __fadd_rn(Su11.f, Stmp2.f);
		Su13.f = __fsub_rn(Su13.f, Stmp1.f);

		Stmp1.f = Ss.f * Su21.f;
		Stmp2.f = Ss.f * Su23.f;
		Su21.f = Sc.f * Su21.f;
		Su23.f = Sc.f * Su23.f;
		Su21.f = __fadd_rn(Su21.f, Stmp2.f);
		Su23.f = __fsub_rn(Su23.f, Stmp1.f);

		Stmp1.f = Ss.f * Su31.f;
		Stmp2.f = Ss.f * Su33.f;
		Su31.f = Sc.f * Su31.f;
		Su33.f = Sc.f * Su33.f;
		Su31.f = __fadd_rn(Su31.f, Stmp2.f);
		Su33.f = __fsub_rn(Su33.f, Stmp1.f);

		// Third Givens Rotation

		Ssh.f = Sa32.f * Sa32.f;
		Ssh.ui = (Ssh.f >= gsmall_number) ? 0xffffffff : 0;
		Ssh.ui = Ssh.ui & Sa32.ui;

		Stmp5.f = 0.f;
		Sch.f = __fsub_rn(Stmp5.f, Sa22.f);
		Sch.f = max(Sch.f, Sa22.f);
		Sch.f = max(Sch.f, gsmall_number);
		Stmp5.ui = (Sa22.f >= Stmp5.f) ? 0xffffffff : 0;

		Stmp1.f = Sch.f * Sch.f;
		Stmp2.f = Ssh.f * Ssh.f;
		Stmp2.f = __fadd_rn(Stmp1.f, Stmp2.f);
		Stmp1.f = __frsqrt_rn(Stmp2.f);

		Stmp4.f = Stmp1.f * 0.5f;
		Stmp3.f = Stmp1.f * Stmp4.f;
		Stmp3.f = Stmp1.f * Stmp3.f;
		Stmp3.f = Stmp2.f * Stmp3.f;
		Stmp1.f = __fadd_rn(Stmp1.f, Stmp4.f);
		Stmp1.f = __fsub_rn(Stmp1.f, Stmp3.f);
		Stmp1.f = Stmp1.f * Stmp2.f;

		Sch.f = __fadd_rn(Sch.f, Stmp1.f);

		Stmp1.ui = ~Stmp5.ui & Ssh.ui;
		Stmp2.ui = ~Stmp5.ui & Sch.ui;
		Sch.ui = Stmp5.ui & Sch.ui;
		Ssh.ui = Stmp5.ui & Ssh.ui;
		Sch.ui = Sch.ui | Stmp1.ui;
		Ssh.ui = Ssh.ui | Stmp2.ui;

		Stmp1.f = Sch.f * Sch.f;
		Stmp2.f = Ssh.f * Ssh.f;
		Stmp2.f = __fadd_rn(Stmp1.f, Stmp2.f);
		Stmp1.f = __frsqrt_rn(Stmp2.f);

		Stmp4.f = Stmp1.f * 0.5f;
		Stmp3.f = Stmp1.f * Stmp4.f;
		Stmp3.f = Stmp1.f * Stmp3.f;
		Stmp3.f = Stmp2.f * Stmp3.f;
		Stmp1.f = __fadd_rn(Stmp1.f, Stmp4.f);
		Stmp1.f = __fsub_rn(Stmp1.f, Stmp3.f);

		Sch.f = Sch.f * Stmp1.f;
		Ssh.f = Ssh.f * Stmp1.f;

		Sc.f = Sch.f * Sch.f;
		Ss.f = Ssh.f * Ssh.f;
		Sc.f = __fsub_rn(Sc.f, Ss.f);
		Ss.f = Ssh.f * Sch.f;
		Ss.f = __fadd_rn(Ss.f, Ss.f);

		//###########################################################
		// Rotate matrix A
		//###########################################################

		Stmp1.f = Ss.f * Sa21.f;
		Stmp2.f = Ss.f * Sa31.f;
		Sa21.f = Sc.f * Sa21.f;
		Sa31.f = Sc.f * Sa31.f;
		Sa21.f = __fadd_rn(Sa21.f, Stmp2.f);
		Sa31.f = __fsub_rn(Sa31.f, Stmp1.f);

		Stmp1.f = Ss.f * Sa22.f;
		Stmp2.f = Ss.f * Sa32.f;
		Sa22.f = Sc.f * Sa22.f;
		Sa32.f = Sc.f * Sa32.f;
		Sa22.f = __fadd_rn(Sa22.f, Stmp2.f);
		Sa32.f = __fsub_rn(Sa32.f, Stmp1.f);

		Stmp1.f = Ss.f * Sa23.f;
		Stmp2.f = Ss.f * Sa33.f;
		Sa23.f = Sc.f * Sa23.f;
		Sa33.f = Sc.f * Sa33.f;
		Sa23.f = __fadd_rn(Sa23.f, Stmp2.f);
		Sa33.f = __fsub_rn(Sa33.f, Stmp1.f);

		//###########################################################
		// Update matrix U
		//###########################################################

		Stmp1.f = Ss.f * Su12.f;
		Stmp2.f = Ss.f * Su13.f;
		Su12.f = Sc.f * Su12.f;
		Su13.f = Sc.f * Su13.f;
		Su12.f = __fadd_rn(Su12.f, Stmp2.f);
		Su13.f = __fsub_rn(Su13.f, Stmp1.f);

		Stmp1.f = Ss.f * Su22.f;
		Stmp2.f = Ss.f * Su23.f;
		Su22.f = Sc.f * Su22.f;
		Su23.f = Sc.f * Su23.f;
		Su22.f = __fadd_rn(Su22.f, Stmp2.f);
		Su23.f = __fsub_rn(Su23.f, Stmp1.f);

		Stmp1.f = Ss.f * Su32.f;
		Stmp2.f = Ss.f * Su33.f;
		Su32.f = Sc.f * Su32.f;
		Su33.f = Sc.f * Su33.f;
		Su32.f = __fadd_rn(Su32.f, Stmp2.f);
		Su33.f = __fsub_rn(Su33.f, Stmp1.f);

		v11 = Sv11.f;
		v12 = Sv12.f;
		v13 = Sv13.f;
		v21 = Sv21.f;
		v22 = Sv22.f;
		v23 = Sv23.f;
		v31 = Sv31.f;
		v32 = Sv32.f;
		v33 = Sv33.f;

		u11 = Su11.f;
		u12 = Su12.f;
		u13 = Su13.f;
		u21 = Su21.f;
		u22 = Su22.f;
		u23 = Su23.f;
		u31 = Su31.f;
		u32 = Su32.f;
		u33 = Su33.f;

		s11 = Sa11.f;
		//s12 = Sa12.f; s13 = Sa13.f; s21 = Sa21.f;
		s22 = Sa22.f;
		//s23 = Sa23.f; s31 = Sa31.f; s32 = Sa32.f;
		s33 = Sa33.f;
	}
}

#endif //_GPU_EIGEN_LIB_CUH_