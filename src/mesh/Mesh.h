#pragma once
#include "Eigen/Eigen"
#include "collision/ConstraintID.h"
#include <vector>
#include <array>

struct mesh_circle {
	double vertex[9][2] = { {0,0},
		{1,0},
	{0.7071067811865476, 0.7071067811865476},
	{0, 1},
	{-0.7071067811865476, 0.7071067811865476},
	{-1, 0},
	{-0.7071067811865476, -0.7071067811865476},
	{0, -1},
	{0.7071067811865476, -0.7071067811865476} };

	//double force[9][2] = {0};
	//double velocity[9][2] = {0};

	int index[8][3] = { {1, 2, 0},
	{2, 3, 0},
	{3, 4, 0},
	{4, 5, 0},
	{5, 6, 0},
	{6, 7, 0},
	{7, 8, 0},
	{8, 1, 0} };
	//Eigen::Matrix2d DM_tetrahedra_inverse[8];
	int vertexNum = 9;
	int triangleNum = 8;
};

struct mesh_rectangle {
	double vertex[15][2] = { {0,1},
		{1,1},
		{1, 2},
		{0, 2},
		{-1, 2},
		{-1, 1},
		{-1, 0},
		{0, 0},
		{1, 0},
		{-1,-1},
		{0,-1},
		{1,-1},
		{-1,-2},
		{0, -2},
		{1,-2} };

	//double force[15][2] = { 0 };
	//double velocity[15][2] = { 0 };

	int index[16][3] = { {1, 2, 0},
		{2, 3, 0},
		{3, 4, 0},
		{4, 5, 0},
		{5, 6, 0},
		{6, 7, 0},
		{7, 8, 0},
		{8, 1, 0},
		{6, 7, 10},
		{7, 8, 10},
		{8, 11, 10},
		{11, 14, 10},
		{13, 14, 10},
		{12, 13, 10},
		{9, 12, 10},
		{6, 9, 10} };
	//Eigen::Matrix2d DM_tetrahedra_inverse[16];
	int vertexNum = 15;
	int triangleNum = 16;
};

struct mesh_cuboid {
	double vertex[8][3] = {
		{-1,1,1},
		{-1,1,-1},
		{1,1,-1},
		{1,1,1},
		{-1,-1,1},
		{-1,-1,-1},
		{1,-1,-1},
		{1,-1,1} };

	//double force[8][3] = { 0 };
	//double velocity[8][3] = { 0 };

	uint64_t index[5][4] = { {0, 4, 7, 5},
		{2, 3, 7, 0},
		{2, 6, 7, 5},
		{0, 1, 2, 5},
		{0, 2, 7, 5} };
	//Eigen::Matrix3d DM_tetrahedra_inverse[5];
	int vertexNum = 8;
	int tetrahedraNum = 5;
};

struct QuadBendingInfo
{
	std::array<size_t, 4> verts;            // [edge1, edge2, trig1 outer, trig2 outer]
	Eigen::Matrix4d Q;						// precomputed local const Hessian

	QuadBendingInfo(size_t edge1, size_t edge2, size_t trig1_outer, size_t trig2_outer, Eigen::Matrix4d Q)
		: Q(Q)
	{
		verts[0] = edge1;
		verts[1] = edge2;
		verts[2] = trig1_outer;
		verts[3] = trig2_outer;
	};
};

class mesh2D {
public:
	std::vector<double> areas;
	std::vector<double> masses;
	std::vector<Eigen::Vector2d> vertexes;
	std::vector<Eigen::Vector3i> triangles;
	std::vector<Eigen::Vector2d> forces;
	std::vector<Eigen::Vector2d> velocities;
	std::vector<Eigen::Matrix2d> DM_triangle_inverse;
	int vertexNum;
	int triangleNum;
	void InitMesh(int type, double scale);
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
	double bendYoungModulus;
	double stretchStiffness;
	double shearStiffness;
	double strainRate;
	double bendingStiffness;
	double maxVolum;
	std::vector<Eigen::Matrix3d> Constraints;
	//std::vector<bool> isDelete;
	//std::vector<int> idmap;
	//std::vector<int> idinverseMap;
	std::vector<double> volum;
	std::vector<double> areas;
	std::vector<double> masses;
	double meanMass;
	double restSNKE;
	double averageEdgeLength;
	double Hhat;
    double Fhat;
	double drag_coeff;
	double IPC_dt;
	double Kappa;
	double bboxDiagSize2;
	double dTol;
	double Newton_Solver_Threshold;
	bool use_barrier;

	std::vector<Eigen::Vector3d> vertexes;
	std::vector<int> boundary_vertexes_indices;
	std::vector<Eigen::Vector3d> v_rest;
	std::vector<Eigen::Vector4i> tetrahedras;
	std::vector<Eigen::Vector3i> triangles;
	std::vector<Eigen::Vector3d> velocities;
	std::vector<Eigen::Matrix3d> DM_tetrahedra_inverse;
	std::vector<Eigen::Matrix2d> DM_triangle_inverse;
	std::vector<Eigen::Matrix<double, 9, 12>> tetPFPX;
	std::vector<Eigen::Matrix<double, 6, 9>> triPFPX;
	std::vector<Eigen::Vector4i> surface;
	std::vector<uint64_t> surfVerts;
	std::vector<std::pair<uint64_t, uint64_t>> surfEdges;
	std::vector<double> Self_lambda_lastH;
	std::vector<double> Environment_lambda_lastH;
	std::vector<Eigen::Vector2d> MMDistCoord;
	std::vector<Eigen::Matrix<double, 3, 2>> MMTanBasis;
	std::vector<int> Environment_ActiveSet, Environment_activeSet_lastH;
	std::vector<MMCVID> Self_ActiveSet, Self_activeSet_lastH;
	std::vector<MMCVID> Self_EE_ActiveSet;
	std::vector<std::pair<int, int>> Self_EEeIe_ActiveSet;
	std::vector<std::pair<int, int>> Self_CCD_ActiveSet;
	std::vector<int> closeConstraintID;
	std::vector<MMCVID> closeMConstraintID;
	std::vector<double> closeConstraintVal;
	std::vector<double> closeMConstraintVal;
	Eigen::Vector3d minCorner, maxCorner;
	Eigen::Vector3d objMinConer, objMaxConer;
	std::vector<int> boundaryTypes;
    std::vector<int> specialPointsArray;
	//IPC
	std::vector<Eigen::Vector3d> xTilta, dx_Elastic, acceleration;
	std::vector<Eigen::Vector3d> V_prev;

	std::vector<Eigen::VectorXd> EKF;

	std::vector<Eigen::Vector2i> tri_edges_adj_points;
	std::vector<Eigen::Vector2i> tri_edges;

	std::vector<QuadBendingInfo> quadBendingInfo;
	std::function<Eigen::Vector3d(Eigen::Vector3d vertex, int step_id, double ipc_dt)> update_hard_constraint_functor = nullptr;
	bool apply_gravity = true;
	bool is_quasi_static = false;
	int D12x12Num;
	int D9x9Num;
	int D6x6Num;
	int D3x3Num;
	int vertexNum;
	int tetrahedraNum;
	int triangleNum;
	void InitMesh(int type, double scale);
	bool load_tetrahedraMesh(const std::string& filename, double scale, Eigen::Vector3d offset);
	bool load_triangleMesh(const std::string& filename, double scale, Eigen::Vector3d position_offset, int type=0);
	void load_test(double scale, int num = 1);
	void getSurface();
	bool output_tetrahedraMesh(const std::string& filename);
	bool output_tetTempData();
	bool load_tetTempData();

	std::vector<Eigen::VectorXd> get_dXn1_dXn();
	std::vector<Eigen::VectorXd> get_dXn1_dVn();
	std::vector<Eigen::VectorXd> get_dVn1_dXn();
	std::vector<Eigen::VectorXd> get_dVn1_dVn();

};

class mesh_obj {
public:
	std::vector<Eigen::Vector3d> vertexes;
	std::vector<Eigen::Vector3d> normals;
	std::vector<Eigen::Vector3i> facenormals;
	std::vector<Eigen::Vector3i> faces;
	int vertexNum;
	int faceNum;
	void InitMesh(int type, double scale);
	bool load_mesh(const std::string& filename, double scale);
};

class model_obj {
public:
	std::vector<mesh_obj> meshes;
	std::vector<std::string> names;
	bool load_model(const std::string& filename);
};

class model_tet {
public:
	std::vector<mesh3D> mesh3Ds;
	//std::vector<std::vector<Eigen::Vector4i>> surfaces;
	std::vector<std::string> names;
	bool load_model(const std::string& filename, int offset);
	void calculate_surface();
};

class fiber_obj {
public:
	std::vector<mesh_obj> muscles;
	std::vector<mesh_obj> tendonIns;
	std::vector<mesh_obj> tendonOuts;
	//std::vector<std::string> names;
	bool load_model(const std::string* filename);
};
