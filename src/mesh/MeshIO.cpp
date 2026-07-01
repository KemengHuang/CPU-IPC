#include "mesh/MeshIO.h"
#include "mesh/Mesh.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <cstring>

using namespace std;
using namespace Eigen;

namespace {

void split(string str, vector<string>& v, string spacer) {
    int pos1, pos2;
    int len = spacer.length();
    pos1 = 0;
    pos2 = str.find(spacer);
    while (pos2 != string::npos) {
        v.push_back(str.substr(pos1, pos2 - pos1));
        pos1 = pos2 + len;
        pos2 = str.find(spacer, pos1);
    }
    if (pos1 != str.length())
        v.push_back(str.substr(pos1));
}

} // namespace

namespace MeshIO {

bool loadTetrahedraMesh(mesh3D& mesh, const std::string& filename, double scale, Eigen::Vector3d offset) {
    ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "unable to read file " << filename << std::endl;
        ifs.close();
        return false;
    }

    double x, y, z;
    int index0, index1, index2, index3;
    string line = "";
    int nodeNumber = 0;
    int elementNumber = 0;
    int offsetVerts = mesh.vertexes.size();
    Vector3d tempMinCorner, tempMaxCorner;

    while (getline(ifs, line)) {
        if (line.length() <= 1) continue;

        if (line.substr(1, 5) == "Nodes") {
            getline(ifs, line);
            nodeNumber = atoi(line.c_str());
            mesh.vertexNum = nodeNumber;

            double xmin = 1e32, ymin = 1e32, zmin = 1e32;
            double xmax = -1e32, ymax = -1e32, zmax = -1e32;
            for (int i = 0; i < nodeNumber; i++) {
                getline(ifs, line);
                vector<std::string> nodePos;
                std::string spacer = " ";
                split(line, nodePos, spacer);
                int posNum = nodePos.size();
                x = atof(nodePos[posNum - 3].c_str());
                y = atof(nodePos[posNum - 2].c_str());
                z = atof(nodePos[posNum - 1].c_str());
                Vector3d vertex = scale * Vector3d(x, y, z) + offset;
                Matrix3d Constraint;
                Constraint.setIdentity();

                Vector3d velocity = Vector3d(0, 0, 0);

                double mass = 0;
                mesh.vertexes.push_back(vertex);

                mesh.velocities.push_back(velocity);
                mesh.Constraints.push_back(Constraint);

                mesh.masses.push_back(mass);

                int boundaryType = 0;
                mesh.boundaryTypes.push_back(boundaryType);
                Vector3d pos = vertex;
                if (xmin > pos[0]) xmin = pos[0];
                if (ymin > pos[1]) ymin = pos[1];
                if (zmin > pos[2]) zmin = pos[2];
                if (xmax < pos[0]) xmax = pos[0];
                if (ymax < pos[1]) ymax = pos[1];
                if (zmax < pos[2]) zmax = pos[2];
            }
            tempMinCorner = Vector3d(xmin, ymin, zmin);
            tempMaxCorner = Vector3d(xmax, ymax, zmax);

            if (mesh.maxCorner[0] < tempMaxCorner[0]) mesh.maxCorner[0] = tempMaxCorner[0];
            if (mesh.maxCorner[1] < tempMaxCorner[1]) mesh.maxCorner[1] = tempMaxCorner[1];
            if (mesh.maxCorner[2] < tempMaxCorner[2]) mesh.maxCorner[2] = tempMaxCorner[2];
            if (mesh.minCorner[0] > tempMinCorner[0]) mesh.minCorner[0] = tempMinCorner[0];
            if (mesh.minCorner[1] > tempMinCorner[1]) mesh.minCorner[1] = tempMinCorner[1];
            if (mesh.minCorner[2] > tempMinCorner[2]) mesh.minCorner[2] = tempMinCorner[2];
        }

        if (line.substr(1, 8) == "Elements") {
            getline(ifs, line);
            elementNumber = atoi(line.c_str());
            mesh.tetrahedraNum = elementNumber;
            for (int i = 0; i < elementNumber; i++) {
                getline(ifs, line);

                vector<std::string> elementIndexex;
                std::string spacer = " ";
                split(line, elementIndexex, spacer);
                int eleNum = elementIndexex.size();
                index0 = atoi(elementIndexex[eleNum - 4].c_str()) - 1;
                index1 = atoi(elementIndexex[eleNum - 3].c_str()) - 1;
                index2 = atoi(elementIndexex[eleNum - 2].c_str()) - 1;
                index3 = atoi(elementIndexex[eleNum - 1].c_str()) - 1;

                Vector4i tetrahedra;
                tetrahedra[0] = (index0 + offsetVerts);
                tetrahedra[1] = (index1 + offsetVerts);
                tetrahedra[2] = (index2 + offsetVerts);
                tetrahedra[3] = (index3 + offsetVerts);

                mesh.tetrahedras.push_back(tetrahedra);
            }
            break;
        }
    }

    if ((tempMaxCorner - tempMinCorner).norm() > (mesh.objMaxCorner - mesh.objMinCorner).norm()) {
        mesh.objMaxCorner = tempMaxCorner;
        mesh.objMinCorner = tempMinCorner;
    }

    ifs.close();
    mesh.V_prev = mesh.vertexes;

    mesh.D12x12Num = 0;
    mesh.D9x9Num = 0;
    mesh.D6x6Num = 0;
    mesh.D3x3Num = 0;
    return true;
}

bool loadTriangleMesh(mesh3D& mesh, const std::string& filename, double scale, Eigen::Vector3d offset, int type) {
    int offsetVerts = mesh.vertexes.size();
    ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "unable to read file " << filename << std::endl;
        ifs.close();
        return false;
    }
    char buffer[1024];
    string line = "";
    int nodeNumber = 0;
    int elementNumber = 0;
    double x, y, z;

    double xmin = 1e32, ymin = 1e32, zmin = 1e32;
    double xmax = -1e32, ymax = -1e32, zmax = -1e32;
    while (getline(ifs, line)) {
        if (line.length() <= 1) continue;
        string key = line.substr(0, 2);
        stringstream ss(line.substr(2));
        if (key == "v ") {
            ss >> x >> y >> z;
            Vector3d vertex = scale * Vector3d(x, y, z) + offset;
            Matrix3d Constraint;
            Constraint.setIdentity();
            Vector3d velocity = Vector3d(0, 0, 0);
            double mass = 0;
            int boundaryType = 0;
            if (type == 2) {
                Constraint.setZero();
                boundaryType = 2;
            }
            mesh.vertexes.push_back(vertex);
            mesh.velocities.push_back(velocity);
            mesh.Constraints.push_back(Constraint);
            mesh.masses.push_back(mass);

            mesh.boundaryTypes.push_back(boundaryType);
            Vector3d pos = vertex;
            if (xmin > pos[0]) xmin = pos[0];
            if (ymin > pos[1]) ymin = pos[1];
            if (zmin > pos[2]) zmin = pos[2];
            if (xmax < pos[0]) xmax = pos[0];
            if (ymax < pos[1]) ymax = pos[1];
            if (zmax < pos[2]) zmax = pos[2];

        } else if (key == "vn") {
            ss >> x >> y >> z;
            Vector3d normal = Vector3d(x, y, z);
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

                if (sscanf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d", &faceVertIndex[0], &uv0, &faceNormalIndex[0],
                           &faceVertIndex[1], &uv1, &faceNormalIndex[1],
                           &faceVertIndex[2], &uv2, &faceNormalIndex[2]) == 9) {
                    faceVertIndex[0] -= 1;
                    faceVertIndex[1] -= 1;
                    faceVertIndex[2] -= 1;

                    faceVertIndex[0] += offsetVerts;
                    faceVertIndex[1] += offsetVerts;
                    faceVertIndex[2] += offsetVerts;
                } else if (sscanf(buffer, "f %d %d %d", &faceVertIndex[0], &faceVertIndex[1], &faceVertIndex[2]) == 3) {
                    faceVertIndex[0] -= 1;
                    faceVertIndex[1] -= 1;
                    faceVertIndex[2] -= 1;

                    faceVertIndex[0] += offsetVerts;
                    faceVertIndex[1] += offsetVerts;
                    faceVertIndex[2] += offsetVerts;
                }

                if (type == 3) {
                    mesh.surface.push_back(Vector4i(faceVertIndex[0], faceVertIndex[1], faceVertIndex[2], 0));
                } else
                    mesh.triangles.push_back(faceVertIndex);
            }
        }
    }
    mesh.vertexNum = mesh.vertexes.size();
    Vector3d tempMinCorner, tempMaxCorner;
    tempMinCorner = Vector3d(xmin, ymin, zmin);
    tempMaxCorner = Vector3d(xmax, ymax, zmax);
    if (mesh.maxCorner[0] < tempMaxCorner[0]) mesh.maxCorner[0] = tempMaxCorner[0];
    if (mesh.maxCorner[1] < tempMaxCorner[1]) mesh.maxCorner[1] = tempMaxCorner[1];
    if (mesh.maxCorner[2] < tempMaxCorner[2]) mesh.maxCorner[2] = tempMaxCorner[2];
    if (mesh.minCorner[0] > tempMinCorner[0]) mesh.minCorner[0] = tempMinCorner[0];
    if (mesh.minCorner[1] > tempMinCorner[1]) mesh.minCorner[1] = tempMinCorner[1];
    if (mesh.minCorner[2] > tempMinCorner[2]) mesh.minCorner[2] = tempMinCorner[2];

    std::set<std::pair<int, int>> edge_set;
    std::map<std::pair<int, int>, std::vector<int>> edge_map;
    std::vector<Eigen::Vector2i> my_edges;
    for (auto tri : mesh.triangles) {
        auto x = tri[0];
        auto y = tri[1];
        auto z = tri[2];
        if (x < y) {
            edge_set.insert(make_pair(x, y));
            edge_map[make_pair(x, y)].emplace_back(z);
        } else {
            edge_set.insert(make_pair(y, x));
            edge_map[make_pair(y, x)].emplace_back(z);
        }

        if (y < z) {
            edge_set.insert(make_pair(y, z));
            edge_map[make_pair(y, z)].emplace_back(x);
        } else {
            edge_set.insert(make_pair(z, y));
            edge_map[make_pair(z, y)].emplace_back(x);
        }

        if (x < z) {
            edge_set.insert(make_pair(x, z));
            edge_map[make_pair(x, z)].emplace_back(y);
        } else {
            edge_set.insert(make_pair(z, x));
            edge_map[make_pair(z, x)].emplace_back(y);
        }
    }

    std::vector<std::vector<int>> temp_edges_adj_points;
    for (auto p : edge_set) {
        if (edge_map[p].size() != 2) continue;
        my_edges.emplace_back(p.first, p.second);
        temp_edges_adj_points.emplace_back(edge_map[make_pair(p.first, p.second)]);
    }
    for (int i = 0; i < temp_edges_adj_points.size(); i++) {
        mesh.tri_edges.emplace_back(Vector2i(my_edges[i].x(), my_edges[i].y()));
        if (temp_edges_adj_points[i].size() == 2)
            mesh.tri_edges_adj_points.emplace_back(Vector2i(temp_edges_adj_points[i][0], temp_edges_adj_points[i][1]));
        else
            mesh.tri_edges_adj_points.emplace_back(Vector2i(temp_edges_adj_points[i][0], -1));
    }

    return true;
}

bool outputTetrahedraMesh(const mesh3D& mesh, const std::string& filename) {
    std::ofstream outmsh1(filename);
    if (!outmsh1) {
        std::cerr << "unable to write file " << filename << std::endl;
        return false;
    }
    outmsh1 << "$Nodes\n";
    outmsh1 << mesh.vertexNum << endl;
    for (int i = 0; i < mesh.vertexNum; i++) {
        outmsh1 << i + 1 << " " << mesh.vertexes[i][0] << " " <<
            mesh.vertexes[i][1] << " " <<
            mesh.vertexes[i][2] << endl;
    }
    outmsh1 << "$Elements\n";
    outmsh1 << mesh.tetrahedraNum << endl;
    for (int i = 0; i < mesh.tetrahedraNum; i++) {
        outmsh1 << i + 1 << " 4 0 " << mesh.tetrahedras[i][0] << " " <<
            mesh.tetrahedras[i][1] << " " <<
            mesh.tetrahedras[i][2] << " " <<
            mesh.tetrahedras[i][3] << endl;
    }
    outmsh1.close();
    return true;
}

bool loadTetTempData(mesh3D& mesh) {
    std::string fileVertex = "tempData/vertex.txt";
    std::string fileXtileVertex = "tempData/vertexXtile.txt";
    ifstream ifs1(fileVertex);
    ifstream ifs2(fileXtileVertex);

    if (!ifs1 || !ifs2) {
        ifs1.close();
        ifs2.close();
        return false;
    }
    double x, y, z;
    int index = 0;
    while (ifs1 >> x >> y >> z) {
        mesh.vertexes[index] = Vector3d(x, y, z);
        mesh.V_prev[index] = mesh.vertexes[index];
        index++;
    }
    index = 0;
    while (ifs2 >> x >> y >> z) {
        mesh.xTilta[index++] = Vector3d(x, y, z);
    }
    ifs1.close();
    ifs2.close();
    return true;
}

bool outputTetTempData(const mesh3D& mesh) {
    std::string fileVertex = "tempData/vertex.txt";
    std::string vertexXtile = "tempData/vertexXtile.txt";

    std::ofstream outmsh1(fileVertex);
    if (!outmsh1) {
        std::cerr << "unable to write file " << fileVertex << std::endl;
        return false;
    }
    for (int i = 0; i < mesh.vertexNum; i++) {
        outmsh1 << mesh.vertexes[i][0] << " " <<
            mesh.vertexes[i][1] << " " <<
            mesh.vertexes[i][2] << endl;
    }
    outmsh1.close();

    std::ofstream outmsh2(vertexXtile);
    if (!outmsh2) {
        std::cerr << "unable to write file " << vertexXtile << std::endl;
        return false;
    }
    for (int i = 0; i < mesh.vertexNum; i++) {
        outmsh2 << mesh.xTilta[i][0] << " " <<
            mesh.xTilta[i][1] << " " <<
            mesh.xTilta[i][2] << endl;
    }
    outmsh2.close();
    return true;
}

} // namespace MeshIO
