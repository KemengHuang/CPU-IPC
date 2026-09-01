#include "SimulationMesh.h"
#include "RuntimePaths.h"
#include <unordered_map>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <set>
#include <stdexcept>
#include <assert.h>

namespace {

struct UnitCubeTetMesh {
    double vertices[8][3] = {
        { -1, 1, 1 }, { -1, 1, -1 }, { 1, 1, -1 }, { 1, 1, 1 },
        { -1, -1, 1 }, { -1, -1, -1 }, { 1, -1, -1 }, { 1, -1, 1 }
    };
    std::uint64_t tetrahedra[5][4] = {
        { 0, 4, 7, 5 }, { 2, 3, 7, 0 }, { 2, 6, 7, 5 },
        { 0, 1, 2, 5 }, { 0, 2, 7, 5 }
    };
    static constexpr int vertexCount = 8;
    static constexpr int tetrahedronCount = 5;
};

} // namespace

unsigned long long calculate_Triangle_hash(const uint64_t &index0, const uint64_t &index1, const uint64_t &index2,
                                           const uint64_t &length) {
    unsigned long long hashVal = (index2 * length + index1) * length + index0;
    return hashVal;
}

class Triangle {

public:
    uint64_t key[3];

    Triangle(const uint64_t *p_key) {
        key[0] = p_key[0];
        key[1] = p_key[1];
        key[2] = p_key[2];
    }

    Triangle(uint64_t key0, uint64_t key1, uint64_t key2) {
        key[0] = key0;
        key[1] = key1;
        key[2] = key2;
    }

    uint64_t operator[](int i) const {
        assert(0 <= i && i <= 2);
        return key[i];
    }

    bool operator<(const Triangle &right) const {
        if (key[0] < right.key[0]) {
            return true;
        } else if (key[0] == right.key[0]) {
            if (key[1] < right.key[1]) {
                return true;
            } else if (key[1] == right.key[1]) {
                if (key[2] < right.key[2]) {
                    return true;
                }
            }
        }
        return false;
    }

    bool operator==(const Triangle &right) const {
        return key[0] == right[0] && key[1] == right[1] && key[2] == right[2];
    }
};

void split(string str, vector<string> &v, string spacer) {
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

void mesh3D::InitMesh(int type, double scale) {
    UnitCubeTetMesh cuboid;
    for (int i = 0; i < UnitCubeTetMesh::vertexCount; i++) {
        Vector3d vertex = scale * Vector3d(
            cuboid.vertices[i][0], cuboid.vertices[i][1], cuboid.vertices[i][2]);
        Vector3d velocity = Vector3d(0, 0, 0);
        double mass = 0;
        vertexes.push_back(vertex);
        velocities.push_back(velocity);
        boundaryTypes.push_back(0);
        masses.push_back(mass);
    }

    for (int i = 0; i < UnitCubeTetMesh::tetrahedronCount; i++) {
        Vector4i tetrahedra;
        for (int j = 0; j < 4; j++) {
            tetrahedra[j] = static_cast<int>(cuboid.tetrahedra[i][j]);
        }
        tetrahedras.push_back(tetrahedra);
    }

    tetrahedraNum = UnitCubeTetMesh::tetrahedronCount;
    vertexNum = UnitCubeTetMesh::vertexCount;
}

void mesh3D::load_test(double scale, int num) {
    if (num < 1 || num > 2) {
        throw std::invalid_argument("load_test supports one or two cuboids");
    }
    vertexNum = 8 * num;
    tetrahedraNum = 5 * num;
    UnitCubeTetMesh test;
    double xmin = 1e32, ymin = 1e32, zmin = 1e32;
    double xmax = -1e32, ymax = -1e32, zmax = -1e32;

    for (int i = 0; i < 5; i++) {
        Vector4i tet;
        for (int j = 0; j < 4; j++) {
            tet[j] = static_cast<int>(test.tetrahedra[i][j]);
        }
        tetrahedras.push_back(tet);
    }

    for (int i = 0; i < 8; i++) {
        Vector3d vert =
                scale * Vector3d(test.vertices[i][0], test.vertices[i][1], test.vertices[i][2]);
        vertexes.push_back(vert);

        Vector3d velocity = Vector3d(0, 0, 0);
        double mass = 0;
        velocities.push_back(velocity);
        boundaryTypes.push_back(0);
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
                    scale * Vector3d(test.vertices[i][0], test.vertices[i][1], test.vertices[i][2])
                    - Vector3d(0.0, 0.5, 0);
            vertexes.push_back(vert);

            Vector3d velocity = Vector3d(0, 0, 0);
            double mass = 0;
            velocities.push_back(velocity);
            boundaryTypes.push_back(0);
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
                tet[j] = static_cast<int>(test.tetrahedra[i][j] + 8);
            }
            tetrahedras.push_back(tet);
        }
    }
    minConer = Vector3d(xmin, ymin, zmin);
    maxConer = Vector3d(xmax, ymax, zmax);
    V_prev = vertexes;
}

bool mesh3D::load_tetrahedraMesh(const std::string &filename, double scale, Vector3d position_offset) {

    ifstream ifs(filename);
    if (!ifs) {
        fprintf(stderr, "unable to read file %s\n", filename.c_str());
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

        if (line.substr(1, 5) == "Nodes") {
            getline(ifs, line);
            nodeNumber = atoi(line.c_str());
            vertexNum = nodeNumber;

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
                Vector3d vertex = scale * Vector3d(x, y, z) + position_offset;
                Vector3d velocity = Vector3d(0, 0, 0);

                double mass = 0;
                vertexes.push_back(vertex);

                velocities.push_back(velocity);

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

        }

        if (line.substr(1, 8) == "Elements") {
            getline(ifs, line);
            elementNumber = atoi(line.c_str());
            tetrahedraNum = elementNumber;
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
                tetrahedra[0] = (index0 + offset);
                tetrahedra[1] = (index1 + offset);
                tetrahedra[2] = (index2 + offset);
                tetrahedra[3] = (index3 + offset);

                tetrahedras.push_back(tetrahedra);
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

    return true;
}

void mesh3D::getSurface() {
    uint64_t length = vertexNum;
    auto triangle_hash = [&](const Triangle &tri) {
        return length * (length * tri[0] + tri[1]) + tri[2];
    };
    //vector<Vector4i> surface;
    std::unordered_map<Triangle, uint64_t, decltype(triangle_hash)> tri2Tet(4 * tetrahedraNum, triangle_hash);
    for (int i = 0; i < tetrahedraNum; i++) {
        const auto &triI = tetrahedras[i];
        for (int j = 0; j < 4; j++) {
            const Triangle &triVInd = Triangle(triI[j % 4], triI[(1 + j) % 4], triI[(2 + j) % 4]);
            if (tri2Tet.find(Triangle(triVInd[0], triVInd[1], triVInd[2])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[0], triVInd[1], triVInd[2])] = tetrahedraNum + 1;
            } else if (tri2Tet.find(Triangle(triVInd[0], triVInd[2], triVInd[1])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[0], triVInd[2], triVInd[1])] = tetrahedraNum + 1;
            } else if (tri2Tet.find(Triangle(triVInd[1], triVInd[0], triVInd[2])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[1], triVInd[0], triVInd[2])] = tetrahedraNum + 1;
            } else if (tri2Tet.find(Triangle(triVInd[1], triVInd[2], triVInd[0])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[1], triVInd[2], triVInd[0])] = tetrahedraNum + 1;
            } else if (tri2Tet.find(Triangle(triVInd[2], triVInd[0], triVInd[1])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[2], triVInd[0], triVInd[1])] = tetrahedraNum + 1;
            } else if (tri2Tet.find(Triangle(triVInd[2], triVInd[1], triVInd[0])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[2], triVInd[1], triVInd[0])] = tetrahedraNum + 1;
            } else {
                tri2Tet[Triangle(triVInd[0], triVInd[1], triVInd[2])] = i;
            }
        }
    }

    for (const auto &triI: tri2Tet) {
        const uint64_t &tetId = triI.second;
        const Triangle &triVInd = triI.first;
        if (tetId < tetrahedraNum) {
            Vector3d vec1 = vertexes[triVInd[1]] - vertexes[triVInd[0]];
            Vector3d vec2 = vertexes[triVInd[2]] - vertexes[triVInd[0]];
            int id3 = 0;
            for (int i = 0; i < 4; i++) {
                if (tetrahedras[tetId][i] != triVInd[0]
                    && tetrahedras[tetId][i] != triVInd[1]
                    && tetrahedras[tetId][i] != triVInd[2]) {
                    id3 = tetrahedras[tetId][i];
                    break;
                }
            }

            Vector3d vec3 = vertexes[id3] - vertexes[triVInd[0]];
            Vector3d n = vec1.cross(vec2);
            if (n.dot(vec3) < 0) {
                surface.push_back(Vector4i(triVInd[0], triVInd[1], triVInd[2], tetId));
            } else {
                surface.push_back(Vector4i(triVInd[0], triVInd[2], triVInd[1], tetId));
            }
        }
    }
    for(const auto & tri:triangles){
        surface.emplace_back(Vector4i(tri[0], tri[1], tri[2], 0));
    }
    vector<bool> flag(vertexNum, false);
    for (const auto &cTri: surface) {
        for (int i = 0; i < 3; i++) {
            if (!flag[cTri[i]]) {
                surfVerts.push_back(cTri[i]);
                flag[cTri[i]] = true;
            }
        }
    }

    set<pair<uint64_t, uint64_t>> SFEdges_set;
    for (const auto &cTri: surface) {
        for (int i = 0; i < 3; i++) {
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

bool mesh3D::output_tetrahedraMesh(const std::string &filename) {
    std::ofstream outmsh1(filename);
    outmsh1 << "$Nodes\n";
    outmsh1 << vertexNum << endl;
    for (int i = 0; i < vertexNum; i++) {
        outmsh1 << i + 1 << " " << vertexes[i][0] << " " <<
                vertexes[i][1] << " " <<
                vertexes[i][2] << endl;
    }
    outmsh1 << "$EndNodes\n";
    outmsh1 << "$Elements\n";
    outmsh1 << tetrahedraNum << endl;
    for (int i = 0; i < tetrahedraNum; i++) {
        outmsh1 << i + 1 << " 4 0 " << tetrahedras[i][0] + 1 << " " <<
                tetrahedras[i][1] + 1 << " " <<
                tetrahedras[i][2] + 1 << " " <<
                tetrahedras[i][3] + 1 << endl;
    }
    outmsh1 << "$EndElements\n";
    return outmsh1.good();
}

bool mesh3D::load_tetTempData() {
    const std::string fileVertex = RuntimePaths::tempFile("vertex.txt");
    const std::string fileXtileVertex = RuntimePaths::tempFile("vertexXtile.txt");
    ifstream ifs1(fileVertex);
    ifstream ifs2(fileXtileVertex);

    if (!ifs1 || !ifs2) {
        return false;
    }

    if (vertexNum < 0
        || vertexes.size() != static_cast<size_t>(vertexNum)
        || V_prev.size() != static_cast<size_t>(vertexNum)
        || inertialTarget.size() != static_cast<size_t>(vertexNum)) {
        return false;
    }

    vector<Vector3d> restoredVertexes(static_cast<size_t>(vertexNum));
    vector<Vector3d> restoredInertialTarget(static_cast<size_t>(vertexNum));
    double x, y, z;

    for (int index = 0; index < vertexNum; ++index) {
        if (!(ifs1 >> x >> y >> z)) {
            return false;
        }
        restoredVertexes[index] = Vector3d(x, y, z);
    }

    for (int index = 0; index < vertexNum; ++index) {
        if (!(ifs2 >> x >> y >> z)) {
            return false;
        }
        restoredInertialTarget[index] = Vector3d(x, y, z);
    }

    // Reject stale state generated for a mesh with a different vertex count.
    if (ifs1 >> x || ifs2 >> x) {
        return false;
    }

    vertexes = restoredVertexes;
    V_prev = restoredVertexes;
    inertialTarget = restoredInertialTarget;
    return true;
}

bool mesh3D::output_tetTempData() {
    if (!RuntimePaths::initialize()) {
        return false;
    }

    std::ofstream outmsh1(RuntimePaths::tempFile("vertex.txt"));
    std::ofstream outmsh2(RuntimePaths::tempFile("vertexXtile.txt"));
    if (!outmsh1 || !outmsh2) {
        return false;
    }

    outmsh1 << std::setprecision(17);
    outmsh2 << std::setprecision(17);
    for (int i = 0; i < vertexNum; i++) {
        outmsh1 << vertexes[i][0] << " " <<
                vertexes[i][1] << " " <<
                vertexes[i][2] << endl;
    }

    for (int i = 0; i < vertexNum; i++) {
        outmsh2 << inertialTarget[i][0] << " " <<
                inertialTarget[i][1] << " " <<
                inertialTarget[i][2] << endl;
    }
    outmsh1.flush();
    outmsh2.flush();
    return outmsh1.good() && outmsh2.good();
}

bool mesh3D::load_triangleMesh(const string &filename, double scale, Vector3d position_offset, int type) {
    int offset = vertexes.size();
    ifstream ifs(filename);
    if (!ifs) {
        fprintf(stderr, "unable to read file %s\n", filename.c_str());
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
            Vector3d vertex = scale * Vector3d(x, y, z) + position_offset;
            Vector3d velocity = Vector3d(0, 0, 0);
            double mass = 0;
            int boundaryType = 0;
            if(type==2){
                boundaryType = 2;
            }
            vertexes.push_back(vertex);
            velocities.push_back(velocity);
            masses.push_back(mass);

            boundaryTypes.push_back(boundaryType);
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
//            normals.push_back(normal);
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

            for (const auto &f: fs) {
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

                    faceVertIndex[0] += offset;
                    faceVertIndex[1] += offset;
                    faceVertIndex[2] += offset;
                } else if (sscanf(buffer, "f %d %d %d", &faceVertIndex[0], &faceVertIndex[1], &faceVertIndex[2]) == 3) {
                    faceVertIndex[0] -= 1;
                    faceVertIndex[1] -= 1;
                    faceVertIndex[2] -= 1;

                    faceVertIndex[0] += offset;
                    faceVertIndex[1] += offset;
                    faceVertIndex[2] += offset;
                }

                if (type == 3) {
                    surface.push_back(Vector4i(faceVertIndex[0],faceVertIndex[1],faceVertIndex[2],0));
                } else
                    triangles.push_back(faceVertIndex);

            }
        }
    }
    vertexNum = vertexes.size();
    Vector3d tempMinConer, tempMaxConer;
    tempMinConer = Vector3d(xmin, ymin, zmin);
    tempMaxConer = Vector3d(xmax, ymax, zmax);
    if (maxConer[0] < tempMaxConer[0]) maxConer[0] = tempMaxConer[0];
    if (maxConer[1] < tempMaxConer[1]) maxConer[1] = tempMaxConer[1];
    if (maxConer[2] < tempMaxConer[2]) maxConer[2] = tempMaxConer[2];
    if (minConer[0] > tempMinConer[0]) minConer[0] = tempMinConer[0];
    if (minConer[1] > tempMinConer[1]) minConer[1] = tempMinConer[1];
    if (minConer[2] > tempMinConer[2]) minConer[2] = tempMinConer[2];









    std::set<std::pair<int, int>> edge_set;
    std::map<std::pair<int, int>, std::vector<int>> edge_map;
    std::vector<Eigen::Vector2i> my_edges;
    for (auto tri : triangles) {
        auto x = tri[0];
        auto y = tri[1];
        auto z = tri[2];
        if (x < y) {
            edge_set.insert(make_pair(x, y));
            edge_map[make_pair(x, y)].emplace_back(z);
        }
        else {
            edge_set.insert(make_pair(y, x));
            edge_map[make_pair(y, x)].emplace_back(z);
        }

        if (y < z) {
            edge_set.insert(make_pair(y, z));
            edge_map[make_pair(y, z)].emplace_back(x);
        }
        else {
            edge_set.insert(make_pair(z, y));
            edge_map[make_pair(z, y)].emplace_back(x);
        }

        if (x < z) {
            edge_set.insert(make_pair(x, z));
            edge_map[make_pair(x, z)].emplace_back(y);
        }
        else {
            edge_set.insert(make_pair(z, x));
            edge_map[make_pair(z, x)].emplace_back(y);
        }
    }

    std::vector<std::vector<int>> temp_edges_adj_points;
    for (auto p : edge_set) {
        if (edge_map[p].size() != 2)continue;
        my_edges.emplace_back(p.first, p.second);
        temp_edges_adj_points.emplace_back(edge_map[make_pair(p.first, p.second)]);
    }
    for (int i = 0; i < temp_edges_adj_points.size(); i++) {
        tri_edges.emplace_back(Vector2i(my_edges[i].x(), my_edges[i].y()));
        if (temp_edges_adj_points[i].size() == 2)
            tri_edges_adj_points.emplace_back(Vector2i(temp_edges_adj_points[i][0], temp_edges_adj_points[i][1]));
        else
            tri_edges_adj_points.emplace_back(Vector2i(temp_edges_adj_points[i][0], -1));
    }

    return true;
}

void SimulationModel::calculateSurface() {

    for (auto &tetMesh: meshes) {
        tetMesh.getSurface();
    }
}
