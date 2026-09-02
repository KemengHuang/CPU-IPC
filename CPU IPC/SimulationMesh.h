#pragma once
#ifndef CIPC_SIMULATION_MESH_H
#define CIPC_SIMULATION_MESH_H
#include "BoundaryConditions.h"
#include "Eigen/Eigen"
#include "EncodedContact.h"
#include "HingeBending.h"
#include <vector>
#include <array>
using namespace std;
using namespace Eigen;

using TetPFPX = Eigen::Matrix<double, 9, 12>;
using TrianglePFPX = Eigen::Matrix<double, 6, 9>;
using TetHessian = Eigen::Matrix<double, 12, 12>;
using TriangleHessian = Eigen::Matrix<double, 9, 9>;
struct QuadBendingInfo
{
	std::array<size_t, 4> verts;            // [edge1, edge2, trig1 outer, trig2 outer]
	Matrix4d Q;								// precomputed local const Hessian
	TetHessian hessianBase;					// Q tensor I3, without material/dt scale

	QuadBendingInfo(size_t edge1, size_t edge2, size_t trig1_outer, size_t trig2_outer, Matrix4d Q)
		: Q(Q), hessianBase(TetHessian::Zero())
	{
		verts[0] = edge1;
		verts[1] = edge2;
		verts[2] = trig1_outer;
		verts[3] = trig2_outer;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				for (int dimension = 0; dimension < 3; ++dimension)
					hessianBase(3 * i + dimension, 3 * j + dimension) = Q(i, j);
	};
};

class mesh3D {
public:
	double density;
	double lengthRateLame;
	double volumeRateLame;
	double lengthRate;
	double volumeRate;
	double friction;

	double cloth_density;
	double YoungModulus;
	double PoissonRate;

	double clothThicness;
	double clothYoungModulus;
	double shearYoungModulus;
	double bendYoungModulus;
	double stretchStiffness;
	double shearStiffness;
	double strainRate;
	double plateRigidity;
	//vector<bool> isDelete;
	//vector<int> idmap;
	//vector<int> idinverseMap;
	vector<double> volum;
	vector<double> areas;
	vector<double> masses;
	double meanMass;
	double averageEdgeLenth;
	double Hhat;
    double Fhat;
	double drag_coeff;
	double IPC_dt;
	double Kappa;
	double bboxDiagSize2;
	double dTol;
	double Newton_Solver_Threshold;
	bool use_barrier;

	vector<Vector3d> vertexes;
	vector<Vector3d> v_rest;
	vector<Vector4i> tetrahedras;
	vector<Vector3i> triangles;
	vector<Vector3d> velocities;
	vector<Matrix3d> DM_tetrahedra_inverse;
	vector<Matrix2d> DM_triangle_inverse;
	vector<TetPFPX> tetrahedraPFPX;
	vector<TrianglePFPX> trianglePFPX;
	vector<Vector4i> surface;
	vector<uint64_t> surfVerts;
	vector<pair<uint64_t, uint64_t>> surfEdges;
	vector<double> Self_lambda_lastH;
	vector<double> Environment_lambda_lastH;
	vector<Eigen::Vector2d> MMDistCoord;
	vector<Eigen::Matrix<double, 3, 2>> MMTanBasis;
	vector<int> Environment_ActiveSet, Environment_activeSet_lastH;
	vector<EncodedContact> Self_ActiveSet, Self_activeSet_lastH;
	vector<EncodedContact> Self_EE_ActiveSet;
	vector<pair<int, int>> Self_EEeIe_ActiveSet;
	vector<pair<int, int>> Self_CCD_ActiveSet;
	vector<int> closeConstraintID;
	std::vector<EncodedContact> closeMConstraintID;
	vector<double> closeConstraintVal;
	vector<double> closeMConstraintVal;
	Vector3d minConer, maxConer;
	Vector3d objMinConer, objMaxConer;
	vector<int> boundaryTypes;
	BoundaryConditionSet boundaryConditions;
	//IPC
	vector<Vector3d> inertialTarget;
	vector<Vector3d> V_prev;

	std::vector<Vector2i> tri_edges_adj_points;
	std::vector<Vector2i> tri_edges;

	vector<QuadBendingInfo> quadBendingInfo;
	vector<HingeBendingInfo> hingeBendingInfo;
	bool apply_gravity = true;
	bool is_quasi_static = false;
	bool resumedFromCheckpoint = false;
	int vertexNum;
	int tetrahedraNum;
	int triangleNum;
	void InitMesh(int type, double scale);
	bool load_tetrahedraMesh(const std::string& filename, double scale, Vector3d offset);
	bool load_triangleMesh(const std::string& filename, double scale, Vector3d position_offset, int type=0);
	void load_test(double scale, int num = 1);
	void getSurface();
	bool output_tetrahedraMesh(const std::string& filename);
	bool output_tetTempData();
	bool load_tetTempData();

};

class SimulationModel {
public:
	vector<mesh3D> meshes;
	void calculateSurface();
};

#endif // CIPC_SIMULATION_MESH_H
