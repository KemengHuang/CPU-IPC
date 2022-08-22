#include "mesh.h"
#include <unordered_map>
#include <fstream>
#include <set>

unsigned long long calculate_Triangle_hash(const uint64_t& index0, const uint64_t& index1, const uint64_t& index2, const uint64_t& length) {
	unsigned long long hashVal = (index2 * length + index1) * length + index0;
	return hashVal;
}

class Triangle {

public:
	uint64_t key[3];

	Triangle(const uint64_t* p_key)
	{
		key[0] = p_key[0];
		key[1] = p_key[1];
		key[2] = p_key[2];
	}
	Triangle(uint64_t key0, uint64_t key1, uint64_t key2)
	{
		key[0] = key0;
		key[1] = key1;
		key[2] = key2;
	}

	uint64_t operator[](int i) const
	{
		assert(0 <= i && i <= 2);
		return key[i];
	}

	bool operator<(const Triangle& right) const
	{
		if (key[0] < right.key[0]) {
			return true;
		}
		else if (key[0] == right.key[0]) {
			if (key[1] < right.key[1]) {
				return true;
			}
			else if (key[1] == right.key[1]) {
				if (key[2] < right.key[2]) {
					return true;
				}
			}
		}
		return false;
	}

	bool operator==(const Triangle& right) const
	{
		return key[0] == right[0] && key[1] == right[1] && key[2] == right[2];
	}
};

void split(string str, vector<string>& v, string spacer)
{
	int pos1, pos2;
	int len = spacer.length();     
	pos1 = 0;
	pos2 = str.find(spacer);
	while (pos2 != string::npos)
	{
		v.push_back(str.substr(pos1, pos2 - pos1));
		pos1 = pos2 + len;
		pos2 = str.find(spacer, pos1);    
	}
	if (pos1 != str.length()) 
		v.push_back(str.substr(pos1));
}

void mesh2D::InitMesh(int type, double scale) {
	if (type == 0) {
		mesh_circle circle;
		for (int i = 0; i < circle.vertexNum; i++) {
			Vector2d vertex = scale * Vector2d(circle.vertex[i][0], circle.vertex[i][1]);
			Vector2d force = Vector2d(0, 0);
			Vector2d velocity = Vector2d(0, 0);
			double mass = 0;
			vertexes.push_back(vertex);
			forces.push_back(force);
			velocities.push_back(velocity);
			masses.push_back(mass);
		}

		for (int i = 0; i < circle.triangleNum; i++) {
			Vector3i triangle;
			triangle[0]=(circle.index[i][0]);
			triangle[1]=(circle.index[i][1]);
			triangle[2]=(circle.index[i][2]);
			//Vector3i triangle = Vector3i(circle.index[i][0], circle.index[i][1], circle.index[i][2]);
			triangles.push_back(triangle);
		}

		triangleNum = circle.triangleNum;
		vertexNum = circle.vertexNum;
	}
	else if (type == 1) {
		mesh_rectangle rectangle;
		for (int i = 0; i < rectangle.vertexNum; i++) {
			Vector2d vertex = scale * Vector2d(rectangle.vertex[i][0], rectangle.vertex[i][1]);
			Vector2d force = Vector2d(0, 0);
			Vector2d velocity = Vector2d(0, 0);
			double mass = 0;
			vertexes.push_back(vertex);
			forces.push_back(force);
			velocities.push_back(velocity);
			masses.push_back(mass);
		}

		for (int i = 0; i < rectangle.triangleNum; i++) {
			Vector3i triangle;
			triangle[0] = (rectangle.index[i][0]);
			triangle[1] = (rectangle.index[i][1]);
			triangle[2] = (rectangle.index[i][2]);
			//Vector3i triangle = Vector3i(circle.index[i][0], circle.index[i][1], circle.index[i][2]);
			triangles.push_back(triangle);
		}

		triangleNum = rectangle.triangleNum;
		vertexNum = rectangle.vertexNum;
	}
}

void mesh3D::InitMesh(int type, double scale) {
	mesh_cuboid cuboid;
	for (int i = 0; i < cuboid.vertexNum; i++) {
		Vector3d vertex = scale * Vector3d(cuboid.vertex[i][0], cuboid.vertex[i][1], cuboid.vertex[i][2]);
		Vector3d force = Vector3d(0, 0, 0);
		Vector3d velocity = Vector3d(0, 0, 0);
		Vector3d d_velocity = Vector3d(0, 0, 0);
		Vector3d d_pos = Vector3d(0, 0, 0);
		double mass = 0;
		Matrix3d Constraint; Constraint.setIdentity();
		vertexes.push_back(vertex);
		//forces.push_back(force);
		velocities.push_back(velocity);
		//d_velocities.push_back(d_velocity);
		Constraints.push_back(Constraint);
		masses.push_back(mass);
		//d_positions.push_back(d_pos);
	}

	for (int i = 0; i < cuboid.tetrahedraNum; i++) {
		Vector4i tetrahedra;
		for (int j = 0; j < 4; j++) {
			tetrahedra[j] = (cuboid.index[i][0]);
		}
		tetrahedras.push_back(tetrahedra);
	}

	tetrahedraNum = cuboid.tetrahedraNum;
	vertexNum = cuboid.vertexNum;
}

void mesh3D::load_test(double scale, int num) {
	vertexNum = 8 * num;
	tetrahedraNum = 5 * num;
	mesh_cuboid test;
	double xmin = DBL_MAX, ymin = DBL_MAX, zmin = DBL_MAX;
	double xmax = -DBL_MAX, ymax = -DBL_MAX, zmax = -DBL_MAX;

	for (int i = 0;i < 5;i++) {
		Vector4i tet;
		for (int j = 0; j < 4; j++) {
			tet[j] = (test.index[i][j] + 8);
		}
		tetrahedras.push_back(tet);

		//tetra_fiberDir.push_back(Vector3d(0, 1, 0));
		//rehabilitate.push_back(-1);
		//isInside.push_back(-1);
		//isJaw.push_back(false);
		//isSkull.push_back(false);
		//isMouth.push_back(false);
		//isMuscle.push_back(false);
	}

	for (int i = 0;i < 8;i++) {
		Vector3d vert = scale*Vector3d(test.vertex[i][0], test.vertex[i][1], test.vertex[i][2]) + Vector3d(0.0, 0.0, 0);
		vertexes.push_back(vert);

		Vector3d d_velocity = Vector3d(0, 0, 0);
		Matrix3d Constraint; Constraint.setIdentity();
		Vector3d force = Vector3d(0, 0, 0);
		Vector3d velocity = Vector3d(0, 0, 0);
		Vector3d d_pos = Vector3d(0, 0, 0);
		double mass = 0;
		//forces.push_back(force);
		velocities.push_back(velocity);
		Constraints.push_back(Constraint);
		//d_velocities.push_back(d_velocity);
		masses.push_back(mass);
		//isDelete.push_back(false);
		//d_positions.push_back(d_pos);
		//externalForce.push_back(Vector3d(0, 0, 0));

		Vector3d pos = vert;
		if (xmin > pos[0]) xmin = pos[0];
		if (ymin > pos[1]) ymin = pos[1];
		if (zmin > pos[2]) zmin = pos[2];
		if (xmax < pos[0]) xmax = pos[0];
		if (ymax < pos[1]) ymax = pos[1];
		if (zmax < pos[2]) zmax = pos[2];
	}
	if (num == 2) {
		for (int i = 0;i < 8;i++) {
			Vector3d vert = scale * Vector3d(test.vertex[i][0], test.vertex[i][1], test.vertex[i][2]) - Vector3d(0.0, 0.5, 0);
			vertexes.push_back(vert);

			Vector3d d_velocity = Vector3d(0, 0, 0);
			Matrix3d Constraint; Constraint.setIdentity();
			Vector3d force = Vector3d(0, 0, 0);
			Vector3d velocity = Vector3d(0, 0, 0);
			Vector3d d_pos = Vector3d(0, 0, 0);
			double mass = 0;
			//forces.push_back(force);
			velocities.push_back(velocity);
			Constraints.push_back(Constraint);
			//d_velocities.push_back(d_velocity);
			masses.push_back(mass);
			//isDelete.push_back(false);
			//d_positions.push_back(d_pos);
			//externalForce.push_back(Vector3d(0, 0, 0));

			Vector3d pos = vert;
			if (xmin > pos[0]) xmin = pos[0];
			if (ymin > pos[1]) ymin = pos[1];
			if (zmin > pos[2]) zmin = pos[2];
			if (xmax < pos[0]) xmax = pos[0];
			if (ymax < pos[1]) ymax = pos[1];
			if (zmax < pos[2]) zmax = pos[2];
		}

		for (int i = 0;i < 5;i++) {
			Vector4i tet;
			for (int j = 0;j < 4;j++) {
				tet[j] = (test.index[i][j] + 8);
			}
			tetrahedras.push_back(tet);

			//tetra_fiberDir.push_back(Vector3d(0, 1, 0));
			//rehabilitate.push_back(-1);
			//isInside.push_back(-1);
			//isJaw.push_back(false);
			//isSkull.push_back(false);
			//isMouth.push_back(false);
			//isMuscle.push_back(false);
		}
	}
	minConer = Vector3d(xmin, ymin, zmin);
	maxConer = Vector3d(xmax, ymax, zmax);
	V_prev = vertexes;
}

bool mesh3D::load_tetrahedraMesh(const std::string& filename, double scale, Vector3d position_offset) {

	ifstream ifs(filename);
	if (!ifs) {

		fprintf(stderr, "unable to read file %s\n", filename.c_str());
		ifs.close();
		exit(-1);
		return false;
	}

	double x, y, z;
	int index0, index1, index2, index3;
	string line = "";
	int nodeNumber = 0;
	int elementNumber = 0;
	int offset = vertexes.size();
	Vector3d tempMinConer, tempMaxConer;
	while (getline(ifs, line)) {
		if (line.length() <= 1) continue;
		if (line == "$Nodes") {
			getline(ifs, line);
			nodeNumber = atoi(line.c_str());
			vertexNum = nodeNumber;

			double xmin = DBL_MAX, ymin = DBL_MAX, zmin = DBL_MAX;
			double xmax = -DBL_MAX, ymax = -DBL_MAX, zmax = -DBL_MAX;
			for (int i = 0; i < nodeNumber; i++) {
				getline(ifs, line);
				vector<std::string> nodePos;
				std::string spacer = " ";
				split(line, nodePos, spacer);
				x = atof(nodePos[1].c_str());
				y = atof(nodePos[2].c_str());
				z = atof(nodePos[3].c_str());
				Vector3d d_velocity = Vector3d(0, 0, 0);
				Vector3d vertex = scale * Vector3d(x, y, z) + position_offset;
				Matrix3d Constraint; Constraint.setIdentity();
				//Vector3d force = Vector3d(0, 0, 0);
				Vector3d velocity = Vector3d(0, 0, 0);
				//Vector3d d_pos = Vector3d(0, 0, 0);
				double mass = 0;
				vertexes.push_back(vertex);
				//forces.push_back(force);
				velocities.push_back(velocity);
				Constraints.push_back(Constraint);
				//d_velocities.push_back(d_velocity);
				masses.push_back(mass);
				//isDelete.push_back(false);
				//d_positions.push_back(d_pos);
				//externalForce.push_back(Vector3d(0, 0, 0));
				int boundaryType = 0;
				boundaryTypes.push_back(boundaryType);
				Vector3d pos = vertex;
				if (xmin > pos[0]) xmin = pos[0];
				if (ymin > pos[1]) ymin = pos[1];
				if (zmin > pos[2]) zmin = pos[2];
				if (xmax < pos[0]) xmax = pos[0];
				if (ymax < pos[1]) ymax = pos[1];
				if (zmax < pos[2]) zmax = pos[2];
			}
			tempMinConer = Vector3d(xmin, ymin, zmin);
			tempMaxConer = Vector3d(xmax, ymax, zmax);

			if (maxConer[0] < tempMaxConer[0]) maxConer[0] = tempMaxConer[0];
			if (maxConer[1] < tempMaxConer[1]) maxConer[1] = tempMaxConer[1];
			if (maxConer[2] < tempMaxConer[2]) maxConer[2] = tempMaxConer[2];
			if (minConer[0] > tempMinConer[0]) minConer[0] = tempMinConer[0];
			if (minConer[1] > tempMinConer[1]) minConer[1] = tempMinConer[1];
			if (minConer[2] > tempMinConer[2]) minConer[2] = tempMinConer[2];

		}

		if (line == "$Elements") {
			getline(ifs, line);
			elementNumber = atoi(line.c_str());
			tetrahedraNum = elementNumber;
			for (int i = 0; i < elementNumber; i++) {
				getline(ifs, line);

				vector<std::string> elementIndexex;
				std::string spacer = " ";
				split(line, elementIndexex, spacer);
				index0 = atoi(elementIndexex[3].c_str()) - 1;
				index1 = atoi(elementIndexex[4].c_str()) - 1;
				index2 = atoi(elementIndexex[5].c_str()) - 1;
				index3 = atoi(elementIndexex[6].c_str()) - 1;

				Vector4i tetrahedra;
				tetrahedra[0] = (index0 + offset);
				tetrahedra[1] = (index1 + offset);
				tetrahedra[2] = (index2 + offset);
				tetrahedra[3] = (index3 + offset);

				tetrahedras.push_back(tetrahedra);

				//tetra_fiberDir.push_back(Vector3d(0, 1, 0));
				//rehabilitate.push_back(-1);
				//isInside.push_back(-1);
				//isJaw.push_back(false);
				//isSkull.push_back(false);
				//isMouth.push_back(false);
				//isMuscle.push_back(false);
			}
			break;
		}
	}

	if ((tempMaxConer - tempMinConer).norm() > (objMaxConer - objMinConer).norm()) {
		objMaxConer = tempMaxConer;
		objMinConer = tempMinConer;
	}

	ifs.close();
	V_prev = vertexes;

	D12x12Num = 0;
	D9x9Num = 0;
	D6x6Num = 0;
	D3x3Num = 0;


	return true;
}



bool mesh3D::load_tetrahedraMesh_IPC_TetMesh(const std::string& filename, Vector3d scale, Vector3d position_offset)
{
	ifstream ifs(filename);
	if (!ifs) {

		fprintf(stderr, "unable to read file %s\n", filename.c_str());
		ifs.close();
		exit(-1);
		return false;
	}

	double x, y, z;
	int index0, index1, index2, index3;
	string line = "";
	int nodeNumber = 0;
	int elementNumber = 0;
	int offset = vertexes.size();
	printf("offset %d\n", offset);
	Vector3d tempMinConer, tempMaxConer;
	while (getline(ifs, line)) {
		if (line.length() <= 1) continue;
		if (line == "$Nodes") {
			getline(ifs, line);
			nodeNumber = atoi(line.c_str());

			double xmin = DBL_MAX, ymin = DBL_MAX, zmin = DBL_MAX;
			double xmax = -DBL_MAX, ymax = -DBL_MAX, zmax = -DBL_MAX;
			for (int i = 0; i < nodeNumber; i++) {
				getline(ifs, line);
				vector<std::string> nodePos;
				std::string spacer = " ";
				split(line, nodePos, spacer);
				x = atof(nodePos[0].c_str());
				y = atof(nodePos[1].c_str());
				z = atof(nodePos[2].c_str());
				Vector3d d_velocity = Vector3d(0, 0, 0);
				Vector3d vertex = Vector3d(scale[0] * x - position_offset.x(), scale[1] * y - position_offset.y(), scale[2] * z - position_offset.z());
				Matrix3d Constraint; Constraint.setIdentity();
				Vector3d force = Vector3d(0, 0, 0);
				Vector3d velocity = Vector3d(0, 0, 0);
				Vector3d d_pos = Vector3d(0, 0, 0);
				double mass = 0;
				vertexes.push_back(vertex);
				//forces.push_back(force);
				velocities.push_back(velocity);
				Constraints.push_back(Constraint);
				//d_velocities.push_back(d_velocity);
				masses.push_back(mass);
				//isDelete.push_back(false);
				//d_positions.push_back(d_pos);
				//externalForce.push_back(Vector3d(0, 0, 0));

				int boundaryType = 0;
				boundaryTypes.push_back(boundaryType);
				Vector3d pos = vertex;
				if (xmin > pos[0]) xmin = pos[0];
				if (ymin > pos[1]) ymin = pos[1];
				if (zmin > pos[2]) zmin = pos[2];
				if (xmax < pos[0]) xmax = pos[0];
				if (ymax < pos[1]) ymax = pos[1];
				if (zmax < pos[2]) zmax = pos[2];
			}
			tempMinConer = Vector3d(xmin, ymin, zmin);
			tempMaxConer = Vector3d(xmax, ymax, zmax);

			if (maxConer[0] < tempMaxConer[0]) maxConer[0] = tempMaxConer[0];
			if (maxConer[1] < tempMaxConer[1]) maxConer[1] = tempMaxConer[1];
			if (maxConer[2] < tempMaxConer[2]) maxConer[2] = tempMaxConer[2];
			if (minConer[0] > tempMinConer[0]) minConer[0] = tempMinConer[0];
			if (minConer[1] > tempMinConer[1]) minConer[1] = tempMinConer[1];
			if (minConer[2] > tempMinConer[2]) minConer[2] = tempMinConer[2];
		}

		if (line == "$Elements") {
			getline(ifs, line);
			elementNumber = atoi(line.c_str());
			for (int i = 0; i < elementNumber; i++) {
				getline(ifs, line);

				vector<std::string> elementIndexex;
				std::string spacer = " ";
				split(line, elementIndexex, spacer);
				index0 = atoi(elementIndexex[1].c_str()) - 1;
				index1 = atoi(elementIndexex[2].c_str()) - 1;
				index2 = atoi(elementIndexex[3].c_str()) - 1;
				index3 = atoi(elementIndexex[4].c_str()) - 1;

				//vector<uint64_t> tetrahedra;
				//tetrahedra.push_back(index0 + offset);
				//tetrahedra.push_back(index1 + offset);
				//tetrahedra.push_back(index2 + offset);
				//tetrahedra.push_back(index3 + offset);
				//tetrahedras.push_back(tetrahedra);

				Vector4i tetrahedra;
				tetrahedra[0] = (index0 + offset);
				tetrahedra[1] = (index1 + offset);
				tetrahedra[2] = (index2 + offset);
				tetrahedra[3] = (index3 + offset);

				tetrahedras.push_back(tetrahedra);
				//tetra_fiberDir.push_back(Vector3d(0, 1, 0));
				//rehabilitate.push_back(-1);
				//isInside.push_back(-1);
				//isJaw.push_back(false);
				//isSkull.push_back(false);
				//isMouth.push_back(false);
				//isMuscle.push_back(false);
			}
			break;
		}
	}
	ifs.close();
	V_prev = vertexes;
	vertexNum = vertexes.size();
	tetrahedraNum = tetrahedras.size();

	if ((tempMaxConer - tempMinConer).norm() > (objMaxConer - objMinConer).norm()) {
		objMaxConer = tempMaxConer;
		objMinConer = tempMinConer;
	}

	D12x12Num = 0;
	D9x9Num = 0;
	D6x6Num = 0;
	D3x3Num = 0;
	return true;

}


bool mesh3D::load_UnitTest(vector<Vector3d> Parray, Vector3d scale, Vector3d position_offset)
{

	double x, y, z;
	int index0, index1, index2, index3;
	int nodeNumber = 0;
	int elementNumber = 0;
	int offset = vertexes.size();

	Vector3d tempMinConer, tempMaxConer;

	//Vector3d P0 = Vector3d(-sqrt(3), 1, -1);
	//Vector3d P1 = Vector3d(sqrt(3), 1, -1);
	//Vector3d P2 = Vector3d(0, 1, 2);
	//Vector3d P3 = Vector3d(0, -1, 0);

	//Vector3d P0 = Vector3d(0, 1, 1);
	//Vector3d P1 = Vector3d(0, 1, -1);
	//Vector3d P2 = Vector3d(-1, -1, 0);
	//Vector3d P3 = Vector3d(1, -1, 0);

	//vector<Vector3d> Parray;
	//Parray.push_back(P0);
	//Parray.push_back(P1);
	//Parray.push_back(P2);
	//Parray.push_back(P3);

	nodeNumber = 4;
	double xmin = DBL_MAX, ymin = DBL_MAX, zmin = DBL_MAX;
	double xmax = -DBL_MAX, ymax = -DBL_MAX, zmax = -DBL_MAX;
	for (auto point : Parray) {

		x = point[0];
		y = point[1];
		z = point[2];
		Vector3d d_velocity = Vector3d(0, 0, 0);
		Vector3d vertex = Vector3d(scale[0] * x - position_offset.x(), scale[1] * y - position_offset.y(), scale[2] * z - position_offset.z());
		Matrix3d Constraint; Constraint.setIdentity();
		Vector3d force = Vector3d(0, 0, 0);
		Vector3d velocity = Vector3d(0, 0, 0);
		Vector3d d_pos = Vector3d(0, 0, 0);
		double mass = 0;
		vertexes.push_back(vertex);
		velocities.push_back(velocity);
		Constraints.push_back(Constraint);
		masses.push_back(mass);

		int boundaryType = 0;
		boundaryTypes.push_back(boundaryType);
		Vector3d pos = vertex;
		if (xmin > pos[0]) xmin = pos[0];
		if (ymin > pos[1]) ymin = pos[1];
		if (zmin > pos[2]) zmin = pos[2];
		if (xmax < pos[0]) xmax = pos[0];
		if (ymax < pos[1]) ymax = pos[1];
		if (zmax < pos[2]) zmax = pos[2];
	}
	tempMinConer = Vector3d(xmin, ymin, zmin);
	tempMaxConer = Vector3d(xmax, ymax, zmax);

	if (maxConer[0] < tempMaxConer[0]) maxConer[0] = tempMaxConer[0];
	if (maxConer[1] < tempMaxConer[1]) maxConer[1] = tempMaxConer[1];
	if (maxConer[2] < tempMaxConer[2]) maxConer[2] = tempMaxConer[2];
	if (minConer[0] > tempMinConer[0]) minConer[0] = tempMinConer[0];
	if (minConer[1] > tempMinConer[1]) minConer[1] = tempMinConer[1];
	if (minConer[2] > tempMinConer[2]) minConer[2] = tempMinConer[2];




	index0 = 0;
	index1 = 1;
	index2 = 2;
	index3 = 3;


	Vector4i tetrahedra;
	tetrahedra[0] = (index0 + offset);
	tetrahedra[1] = (index1 + offset);
	tetrahedra[2] = (index2 + offset);
	tetrahedra[3] = (index3 + offset);

	tetrahedras.push_back(tetrahedra);


	V_prev = vertexes;
	vertexNum = vertexes.size();
	tetrahedraNum = tetrahedras.size();

	if ((tempMaxConer - tempMinConer).norm() > (objMaxConer - objMinConer).norm()) {
		objMaxConer = tempMaxConer;
		objMinConer = tempMinConer;
	}

	D12x12Num = 0;
	D9x9Num = 0;
	D6x6Num = 0;
	D3x3Num = 0;
	return true;

}

bool mesh3D::load_UnitTest2(Vector3d scale, Vector3d position_offset)
{

	double x, y, z;
	int index0, index1, index2, index3;
	int nodeNumber = 0;
	int elementNumber = 0;
	int offset = vertexes.size();

	Vector3d tempMinConer, tempMaxConer;

	Vector3d P0 = Vector3d(-3, 1, 0);
	Vector3d P1 = Vector3d(-3, -1, 0);
	Vector3d P2 = Vector3d(0, 0, -2);
	Vector3d P3 = Vector3d(0, 0, 2);
	Vector3d P4 = Vector3d(3, 1, 0);
	Vector3d P5 = Vector3d(3, -1, 0);


	//Vector3d P0 = Vector3d(0, 1, 1);
	//Vector3d P1 = Vector3d(0, 1, -1);
	//Vector3d P2 = Vector3d(-1, -1, 0);
	//Vector3d P3 = Vector3d(1, -1, 0);

	vector<Vector3d> Parray;
	Parray.push_back(P0);
	Parray.push_back(P1);
	Parray.push_back(P2);
	Parray.push_back(P3);
	Parray.push_back(P4);
	Parray.push_back(P5);

	nodeNumber = 6;
	double xmin = DBL_MAX, ymin = DBL_MAX, zmin = DBL_MAX;
	double xmax = -DBL_MAX, ymax = -DBL_MAX, zmax = -DBL_MAX;
	for (auto point : Parray) {

		x = point[0];
		y = point[1];
		z = point[2];
		Vector3d d_velocity = Vector3d(0, 0, 0);
		Vector3d vertex = Vector3d(scale[0] * x - position_offset.x(), scale[1] * y - position_offset.y(), scale[2] * z - position_offset.z());
		Matrix3d Constraint; Constraint.setIdentity();
		Vector3d force = Vector3d(0, 0, 0);
		Vector3d velocity = Vector3d(0, 0, 0);
		Vector3d d_pos = Vector3d(0, 0, 0);
		double mass = 0;
		vertexes.push_back(vertex);
		velocities.push_back(velocity);
		Constraints.push_back(Constraint);
		masses.push_back(mass);

		int boundaryType = 0;
		boundaryTypes.push_back(boundaryType);
		Vector3d pos = vertex;
		if (xmin > pos[0]) xmin = pos[0];
		if (ymin > pos[1]) ymin = pos[1];
		if (zmin > pos[2]) zmin = pos[2];
		if (xmax < pos[0]) xmax = pos[0];
		if (ymax < pos[1]) ymax = pos[1];
		if (zmax < pos[2]) zmax = pos[2];
	}
	tempMinConer = Vector3d(xmin, ymin, zmin);
	tempMaxConer = Vector3d(xmax, ymax, zmax);

	if (maxConer[0] < tempMaxConer[0]) maxConer[0] = tempMaxConer[0];
	if (maxConer[1] < tempMaxConer[1]) maxConer[1] = tempMaxConer[1];
	if (maxConer[2] < tempMaxConer[2]) maxConer[2] = tempMaxConer[2];
	if (minConer[0] > tempMinConer[0]) minConer[0] = tempMinConer[0];
	if (minConer[1] > tempMinConer[1]) minConer[1] = tempMinConer[1];
	if (minConer[2] > tempMinConer[2]) minConer[2] = tempMinConer[2];




	index0 = 0;
	index1 = 1;
	index2 = 2;
	index3 = 3;


	Vector4i tetrahedra;
	tetrahedra[0] = (index0 + offset);
	tetrahedra[1] = (index1 + offset);
	tetrahedra[2] = (index2 + offset);
	tetrahedra[3] = (index3 + offset);

	tetrahedras.push_back(tetrahedra);



	index0 = 2;
	index1 = 3;
	index2 = 4;
	index3 = 5;


	Vector4i tetrahedra2;
	tetrahedra2[0] = (index0 + offset);
	tetrahedra2[1] = (index1 + offset);
	tetrahedra2[2] = (index2 + offset);
	tetrahedra2[3] = (index3 + offset);

	tetrahedras.push_back(tetrahedra2);


	V_prev = vertexes;
	vertexNum = vertexes.size();
	tetrahedraNum = tetrahedras.size();

	if ((tempMaxConer - tempMinConer).norm() > (objMaxConer - objMinConer).norm()) {
		objMaxConer = tempMaxConer;
		objMinConer = tempMinConer;
	}

	D12x12Num = 0;
	D9x9Num = 0;
	D6x6Num = 0;
	D3x3Num = 0;
	return true;

}


void mesh3D::getSurface() {
	uint64_t length = vertexNum;
	auto triangle_hash = [&](const Triangle& tri) {
		return length * (length * tri[0] + tri[1]) + tri[2];
	};
	//vector<Vector4i> surface;
	std::unordered_map<Triangle, uint64_t, decltype(triangle_hash)> tri2Tet(4 * tetrahedraNum, triangle_hash);
	for (int i = 0;i < tetrahedraNum;i++) {
		const auto& triI = tetrahedras[i];
		for (int j = 0;j < 4;j++) {
			const Triangle& triVInd = Triangle(triI[j % 4], triI[(1 + j) % 4], triI[(2 + j) % 4]);
			if (tri2Tet.find(Triangle(triVInd[0], triVInd[1], triVInd[2])) != tri2Tet.end()) {
				tri2Tet[Triangle(triVInd[0], triVInd[1], triVInd[2])] = tetrahedraNum + 1;
			}
			else if (tri2Tet.find(Triangle(triVInd[0], triVInd[2], triVInd[1])) != tri2Tet.end()) {
				tri2Tet[Triangle(triVInd[0], triVInd[2], triVInd[1])] = tetrahedraNum + 1;
			}
			else if (tri2Tet.find(Triangle(triVInd[1], triVInd[0], triVInd[2])) != tri2Tet.end()) {
				tri2Tet[Triangle(triVInd[1], triVInd[0], triVInd[2])] = tetrahedraNum + 1;
			}
			else if (tri2Tet.find(Triangle(triVInd[1], triVInd[2], triVInd[0])) != tri2Tet.end()) {
				tri2Tet[Triangle(triVInd[1], triVInd[2], triVInd[0])] = tetrahedraNum + 1;
			}
			else if (tri2Tet.find(Triangle(triVInd[2], triVInd[0], triVInd[1])) != tri2Tet.end()) {
				tri2Tet[Triangle(triVInd[2], triVInd[0], triVInd[1])] = tetrahedraNum + 1;
			}
			else if (tri2Tet.find(Triangle(triVInd[2], triVInd[1], triVInd[0])) != tri2Tet.end()) {
				tri2Tet[Triangle(triVInd[2], triVInd[1], triVInd[0])] = tetrahedraNum + 1;
			}
			else {
				tri2Tet[Triangle(triVInd[0], triVInd[1], triVInd[2])] = i;
			}
		}
	}

	for (const auto& triI : tri2Tet) {
		const uint64_t& tetId = triI.second;
		const Triangle& triVInd = triI.first;
		if (tetId < tetrahedraNum) {
			Vector3d vec1 = vertexes[triVInd[1]] - vertexes[triVInd[0]];
			Vector3d vec2 = vertexes[triVInd[2]] - vertexes[triVInd[0]];
			int id3 = 0;
			for (int i = 0;i < 4;i++) {
				if (tetrahedras[tetId][i] != triVInd[0]
					&& tetrahedras[tetId][i] != triVInd[1]
					&& tetrahedras[tetId][i] != triVInd[2]) {
					id3 = tetrahedras[tetId][i]; break;
				}
			}

			Vector3d vec3 = vertexes[id3] - vertexes[triVInd[0]];
			Vector3d n = vec1.cross(vec2);
			if (n.dot(vec3) < 0) {
				surface.push_back(Vector4i(triVInd[0], triVInd[1], triVInd[2], tetId));
			}
			else {
				surface.push_back(Vector4i(triVInd[0], triVInd[2], triVInd[1], tetId));
			}
		}
	}

	vector<bool> flag(vertexNum, false);
	for (const auto& cTri : surface) {
		for (int i = 0;i < 3;i++) {
			if (!flag[cTri[i]]) {
				surfVerts.push_back(cTri[i]);
				flag[cTri[i]] = true;
			}
		}
	}

	set<pair<uint64_t, uint64_t>> SFEdges_set;
	for (const auto& cTri : surface) {
		for (int i = 0;i < 3;i++) {
			if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[1], cTri[0])) == SFEdges_set.end()) {
				SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[0], cTri[1]));
			}
			if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[2], cTri[1])) == SFEdges_set.end()) {
				SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[1], cTri[2]));
			}
			if (SFEdges_set.find(pair<uint64_t, uint64_t>(cTri[0], cTri[2])) == SFEdges_set.end()) {
				SFEdges_set.insert(pair<uint64_t, uint64_t>(cTri[2], cTri[0]));
			}
		}
	}
	surfEdges = vector<pair<uint64_t, uint64_t>>(SFEdges_set.begin(), SFEdges_set.end());
}

bool mesh3D::output_tetrahedraMesh(const std::string& filename) {
	std::ofstream outmsh1(filename);
	outmsh1 << "$Nodes\n";
	outmsh1 << vertexNum << endl;
	for (int i = 0; i < vertexNum; i++) {
		outmsh1 << i + 1 << " " << vertexes[i][0] << " " <<
			vertexes[i][1] << " " <<
			vertexes[i][2] << endl;
	}
	outmsh1 << "$Elements\n";
	outmsh1 << tetrahedraNum << endl;
	for (int i = 0; i < tetrahedraNum; i++) {
		outmsh1 << i + 1 << " 4 0 " << tetrahedras[i][0] << " " <<
			tetrahedras[i][1] << " " <<
			tetrahedras[i][2] << " " <<
			tetrahedras[i][3] << endl;
	}
	outmsh1.close();
	return true;
}

bool mesh3D::load_tetTempData() {
	std::string fileVertex = "tempData/vertex.txt";
	std::string fileXtileVertex = "tempData/vertexXtile.txt";
	ifstream ifs1(fileVertex);
	ifstream ifs2(fileXtileVertex);

	if (!ifs1|| !ifs2) {
		ifs1.close();
		ifs2.close();
		return false;
	}
	double x, y, z;
	int index = 0;
	while (ifs1 >> x >> y >> z) {
		vertexes[index] = Vector3d(x, y, z);
		V_prev[index] = vertexes[index];
		index++;
	}
	index = 0;
	while (ifs2 >> x >> y >> z) {
		xTilta[index++] = Vector3d(x, y, z);
	}
	ifs1.close();
	ifs2.close();
	return true;
}

bool mesh3D::output_tetTempData() {
	std::string fileVertex = "tempData/vertex.txt";
	std::string vertexXtile = "tempData/vertexXtile.txt";


	std::ofstream outmsh1(fileVertex);
	for (int i = 0; i < vertexNum; i++) {
		outmsh1 << vertexes[i][0] << " " <<
			vertexes[i][1] << " " <<
			vertexes[i][2] << endl;
	}
	outmsh1.close();

	std::ofstream outmsh2(vertexXtile);
	for (int i = 0; i < vertexNum; i++) {
		outmsh2 << xTilta[i][0] << " " <<
			xTilta[i][1] << " " <<
			xTilta[i][2] << endl;
	}
	outmsh2.close();
	return true;
}

bool mesh_obj::load_mesh(const std::string& filename, double scale) {
	ifstream ifs(filename);
	if (!ifs) {

		fprintf(stderr, "unable to read file %s\n", filename.c_str());
		ifs.close();
		exit(-1);
		return false;
	}
	char buffer[1024];
	string line = "";
	int nodeNumber = 0;
	int elementNumber = 0;
	double x, y, z;

	while (getline(ifs, line)) {
		string key = line.substr(0, 2);
		stringstream ss(line.substr(2));
		if (key == "v ") {
			ss >> x >> y >> z;
			Vector3d vertex = scale * Vector3d(x, y, z);
			vertexes.push_back(vertex);
		}
		else if (key == "vn") {
			ss >> x >> y >> z;
			Vector3d normal = Vector3d(x, y, z);
			normals.push_back(normal);
		}
		else if (key == "f ") {
			if (line.length() >= 1024) {
				printf("[WARN]: skip line due to exceed max buffer length (1024).\n");
				continue;
			}

			std::vector<string> fs;

			{
				string buf;
				stringstream ss(line);
				vector<string> tokens;
				while (ss >> buf)
					tokens.push_back(buf);

				for (size_t index = 3; index < tokens.size(); index += 1) {
					fs.push_back("f " + tokens[1] + " " + tokens[index - 1] + " " + tokens[index]);
				}
			}

			int uv0, uv1, uv2;

			for (const auto& f : fs) {
				memset(buffer, 0, sizeof(char) * 1024);
				std::copy(f.begin(), f.end(), buffer);

				Vector3i faceVertIndex;
				Vector3i faceNormalIndex;
				
				sscanf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d", &faceVertIndex[0], &uv0, &faceNormalIndex[0],
					&faceVertIndex[1], &uv1, &faceNormalIndex[1],
					&faceVertIndex[2], &uv2, &faceNormalIndex[2]);

				faceVertIndex[0] -= 1;
				faceVertIndex[1] -= 1;
				faceVertIndex[2] -= 1;
				faces.push_back(faceVertIndex);
				facenormals.push_back(faceNormalIndex);

			}
		}
	}

	vertexNum = vertexes.size();
	faceNum = faces.size();

	for (int i = 0; i < faceNum; i++) {
		Vector3d A = vertexes[faces[i][1]] - vertexes[faces[i][0]];
		Vector3d B = vertexes[faces[i][2]] - vertexes[faces[i][0]];
		Vector3d normal = (A.cross(B)).normalized();
		normals.push_back(normal);
	}




	return true;
}

bool model_obj::load_model(const std::string& filename) {
	ifstream ifs(filename);
	if (!ifs) {

		fprintf(stderr, "unable to read file %s\n", filename.c_str());
		ifs.close();
		exit(-1);
		return false;
	}

	string meshName;
	//int id = 0;
	float scale = 1;
	while (ifs >> meshName) {
		mesh_obj mesh;
		mesh.load_mesh(meshName, scale);
		meshes.push_back(mesh);
		names.push_back(meshName);
	}

	return true;
}

bool fiber_obj::load_model(const std::string* filename) {

	ifstream ifsMuscle(filename[0]), ifsTendonIn(filename[1]), ifsTendonOut(filename[2]);
	if (!ifsMuscle|| !ifsTendonIn || !ifsTendonOut) {

		fprintf(stderr, "unable to read file %s\n");
		ifsMuscle.close();
		ifsTendonIn.close(); 
		ifsTendonOut.close();
		exit(-1);
		return false;
	}

	string meshName;
	float scale = 1;
	while (ifsMuscle >> meshName) {
		mesh_obj mesh;
		mesh.load_mesh(meshName, scale);
		muscles.push_back(mesh);
		//names.push_back(meshName);
	}
	string tendonName;
	while (ifsTendonIn >> tendonName) {
		mesh_obj mesh;
		mesh.load_mesh(tendonName, scale);
		tendonIns.push_back(mesh);
		//names.push_back(meshName);
	}

	while (ifsTendonOut >> tendonName) {
		mesh_obj mesh;
		mesh.load_mesh(tendonName, scale);
		tendonOuts.push_back(mesh);
		//names.push_back(meshName);
	}

	return true;
}
#include<iostream>
void model_tet::calculate_surface() {

	for (auto& tetMesh : mesh3Ds) {
		tetMesh.getSurface();
	}
}

bool model_tet::load_model(const std::string& filename, int offset) {
	

	return true;
}