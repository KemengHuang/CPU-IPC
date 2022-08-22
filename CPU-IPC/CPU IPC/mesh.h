#pragma once
#ifndef FEM_MESH_H
#define FEM_MESH_H
#include "Eigen/Eigen"
#include "mIPC.h"
#include <vector>
using namespace std;
using namespace Eigen;
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
	//Matrix2d DM_triangle_inverse[8];
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
	//Matrix2d DM_triangle_inverse[16];
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
	//Matrix3d DM_tetrahedra_inverse[5];
	int vertexNum = 8;
	int tetrahedraNum = 5;
};


class mesh2D {
public:
	vector<double> areas;
	vector<double> masses;
	vector<Vector2d> vertexes;
	vector<Vector3i> triangles;
	vector<Vector2d> forces;
	vector<Vector2d> velocities;
	vector<Matrix2d> DM_triangle_inverse;
	int vertexNum;
	int triangleNum;
	void InitMesh(int type, double scale);
};

class mesh3D {
public:
	double maxVolum;
	vector<Matrix3d> Constraints;
	//vector<bool> isDelete;
	//vector<int> idmap;
	//vector<int> idinverseMap;
	vector<double> volum;
	vector<double> masses;
	double meanMass;
	double restSNKE;
	double averageEdgeLenth;
	double Hhat;
	double Kappa;
	double bboxDiagSize2;
	double dTol;
	vector<Vector3d> vertexes;
	vector<Vector3d> v_rest;
	vector<Vector4i> tetrahedras;
	vector<Vector3d> velocities;
	vector<Matrix3d> DM_triangle_inverse;
	vector<Vector4i> surface;
	vector<uint64_t> surfVerts;
	vector<pair<uint64_t, uint64_t>> surfEdges;
	vector<double> Self_lambda_lastH;
	vector<double> Environment_lambda_lastH;
	vector<Eigen::Vector2d> MMDistCoord;
	vector<Eigen::Matrix<double, 3, 2>> MMTanBasis;
	vector<int> Environment_ActiveSet, Environment_activeSet_lastH;
	vector<MMCVID> Self_ActiveSet, Self_activeSet_lastH;
	vector<MMCVID> Self_EE_ActiveSet;
	vector<pair<int, int>> Self_EEeIe_ActiveSet;
	vector<pair<int, int>> Self_CCD_ActiveSet;
	vector<int> closeConstraintID;
	std::vector<MMCVID> closeMConstraintID;
	vector<double> closeConstraintVal;
	vector<double> closeMConstraintVal;
	Vector3d minConer, maxConer;
	Vector3d objMinConer, objMaxConer;
	vector<int> boundaryTypes;
	//IPC
	vector<Vector3d> xTilta, dx_Elastic, acceleration;
	vector<Vector3d> V_prev;
	int D12x12Num;
	int D9x9Num;
	int D6x6Num;
	int D3x3Num;
	int vertexNum;
	int tetrahedraNum;
	void InitMesh(int type, double scale);
	bool load_tetrahedraMesh(const std::string& filename, double scale, Vector3d offset);
	bool load_tetrahedraMesh_IPC_TetMesh(const std::string& filename, double scale, Vector3d position_offset);
	void load_test(double scale, int num = 1);
	void getSurface();
	bool output_tetrahedraMesh(const std::string& filename);
	bool output_tetTempData();
	bool load_tetTempData();
};

class mesh_obj {
public:
	vector<Vector3d> vertexes;
	vector<Vector3d> normals;
	vector<Vector3i> facenormals;
	vector<Vector3i> faces;
	int vertexNum;
	int faceNum;
	void InitMesh(int type, double scale);
	bool load_mesh(const std::string& filename, double scale);
};

class model_obj {
public:
	vector<mesh_obj> meshes;
	vector<string> names;
	bool load_model(const std::string& filename);
};

class model_tet {
public:
	vector<mesh3D> mesh3Ds;
	//vector<vector<Vector4i>> surfaces;
	vector<string> names;
	bool load_model(const std::string& filename, int offset);
	void calculate_surface();
};

class fiber_obj {
public:
	vector<mesh_obj> muscles;
	vector<mesh_obj> tendonIns;
	vector<mesh_obj> tendonOuts;
	//vector<string> names;
	bool load_model(const std::string* filename);
};

#endif // !FEM_MESH.H