#include "mesh/Mesh.h"
#include "mesh/MeshGeometry.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>

using namespace std;
using namespace Eigen;

void mesh3D::InitMesh(int type, double scale) {
    mesh_cuboid cuboid;
    for (int i = 0; i < cuboid.vertexNum; i++) {
        Vector3d vertex = scale * Vector3d(cuboid.vertex[i][0], cuboid.vertex[i][1], cuboid.vertex[i][2]);
        Vector3d velocity = Vector3d(0, 0, 0);
        Matrix3d Constraint;
        Constraint.setIdentity();
        vertexes.push_back(vertex);
        velocities.push_back(velocity);
        Constraints.push_back(Constraint);
        masses.push_back(0);
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
    double xmin = 1e32, ymin = 1e32, zmin = 1e32;
    double xmax = -1e32, ymax = -1e32, zmax = -1e32;

    for (int i = 0; i < 5; i++) {
        Vector4i tet;
        for (int j = 0; j < 4; j++) {
            tet[j] = (test.index[i][j] + 8);
        }
        tetrahedras.push_back(tet);
    }

    for (int i = 0; i < 8; i++) {
        Vector3d vert =
            scale * Vector3d(test.vertex[i][0], test.vertex[i][1], test.vertex[i][2]) + Vector3d(0.0, 0.0, 0);
        vertexes.push_back(vert);

        Vector3d velocity = Vector3d(0, 0, 0);
        Matrix3d Constraint;
        Constraint.setIdentity();
        double mass = 0;
        velocities.push_back(velocity);
        Constraints.push_back(Constraint);
        masses.push_back(mass);

        Vector3d pos = vert;
        if (xmin > pos[0]) xmin = pos[0];
        if (ymin > pos[1]) ymin = pos[1];
        if (zmin > pos[2]) zmin = pos[2];
        if (xmax < pos[0]) xmax = pos[0];
        if (ymax < pos[1]) ymax = pos[1];
        if (zmax < pos[2]) zmax = pos[2];
    }
    if (num == 2) {
        for (int i = 0; i < 8; i++) {
            Vector3d vert =
                scale * Vector3d(test.vertex[i][0], test.vertex[i][1], test.vertex[i][2]) - Vector3d(0.0, 0.5, 0);
            vertexes.push_back(vert);

            Vector3d velocity = Vector3d(0, 0, 0);
            Matrix3d Constraint;
            Constraint.setIdentity();
            double mass = 0;
            velocities.push_back(velocity);
            Constraints.push_back(Constraint);
            masses.push_back(mass);

            Vector3d pos = vert;
            if (xmin > pos[0]) xmin = pos[0];
            if (ymin > pos[1]) ymin = pos[1];
            if (zmin > pos[2]) zmin = pos[2];
            if (xmax < pos[0]) xmax = pos[0];
            if (ymax < pos[1]) ymax = pos[1];
            if (zmax < pos[2]) zmax = pos[2];
        }

        for (int i = 0; i < 5; i++) {
            Vector4i tet;
            for (int j = 0; j < 4; j++) {
                tet[j] = (test.index[i][j] + 8);
            }
            tetrahedras.push_back(tet);
        }
    }
    minCorner = Vector3d(xmin, ymin, zmin);
    maxCorner = Vector3d(xmax, ymax, zmax);
    V_prev = vertexes;
}

void mesh3D::getSurface() {
    MeshGeometry::computeSurface(*this);
}

bool mesh_obj::load_mesh(const std::string& filename, double scale) {
    ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "unable to read file " << filename << std::endl;
        ifs.close();
        return false;
    }
    char buffer[1024];
    string line = "";
    double x, y, z;

    while (getline(ifs, line)) {
        string key = line.substr(0, 2);
        stringstream ss(line.substr(2));
        if (key == "v ") {
            ss >> x >> y >> z;
            Vector3d vertex = scale * Vector3d(x, y, z);
            vertexes.push_back(vertex);
        } else if (key == "vn") {
            ss >> x >> y >> z;
            Vector3d normal = Vector3d(x, y, z);
            normals.push_back(normal);
        } else if (key == "f ") {
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
        std::cerr << "unable to read file " << filename << std::endl;
        ifs.close();
        return false;
    }

    string meshName;
    float scale = 1;
    while (ifs >> meshName) {
        mesh_obj mesh;
        if (!mesh.load_mesh(meshName, scale)) {
            std::cerr << "failed to load mesh " << meshName << std::endl;
            continue;
        }
        meshes.push_back(mesh);
        names.push_back(meshName);
    }

    return true;
}

void model_tet::calculate_surface() {
    for (auto& tetMesh : mesh3Ds) {
        tetMesh.getSurface();
    }
}

bool model_tet::load_model(const std::string& filename, int offset) {
    return true;
}

vector<Eigen::VectorXd> mesh3D::get_dXn1_dXn() {
    return EKF;
}
vector<Eigen::VectorXd> mesh3D::get_dXn1_dVn() {
    vector<Eigen::VectorXd> dXn1_dVn;
    for (int i = 0; i < vertexNum * 3; i++) {
        dXn1_dVn.push_back(IPC_dt * EKF[i]);
    }
    return dXn1_dVn;
}
vector<Eigen::VectorXd> mesh3D::get_dVn1_dXn() {
    vector<Eigen::VectorXd> dVn1_dXn;
    double one_div_dt = 1.0 / IPC_dt;
    for (int i = 0; i < vertexNum * 3; i++) {
        dVn1_dXn.push_back(one_div_dt * EKF[i]);
    }
    return dVn1_dXn;
}
vector<Eigen::VectorXd> mesh3D::get_dVn1_dVn() {
    return EKF;
}
