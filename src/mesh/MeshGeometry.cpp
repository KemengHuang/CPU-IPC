#include "mesh/MeshGeometry.h"
#include "mesh/Mesh.h"

#include <unordered_map>
#include <set>
#include <vector>

using namespace std;
using namespace Eigen;

namespace MeshGeometry {

void computeSurface(mesh3D& mesh) {
    uint64_t length = mesh.vertexNum;
    auto triangle_hash = [&](const Triangle& tri) {
        return length * (length * tri[0] + tri[1]) + tri[2];
    };
    std::unordered_map<Triangle, uint64_t, decltype(triangle_hash)> tri2Tet(4 * mesh.tetrahedraNum, triangle_hash);
    for (int i = 0; i < mesh.tetrahedraNum; i++) {
        const auto& triI = mesh.tetrahedras[i];
        for (int j = 0; j < 4; j++) {
            const Triangle& triVInd = Triangle(triI[j % 4], triI[(1 + j) % 4], triI[(2 + j) % 4]);
            if (tri2Tet.find(Triangle(triVInd[0], triVInd[1], triVInd[2])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[0], triVInd[1], triVInd[2])] = mesh.tetrahedraNum + 1;
            } else if (tri2Tet.find(Triangle(triVInd[0], triVInd[2], triVInd[1])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[0], triVInd[2], triVInd[1])] = mesh.tetrahedraNum + 1;
            } else if (tri2Tet.find(Triangle(triVInd[1], triVInd[0], triVInd[2])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[1], triVInd[0], triVInd[2])] = mesh.tetrahedraNum + 1;
            } else if (tri2Tet.find(Triangle(triVInd[1], triVInd[2], triVInd[0])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[1], triVInd[2], triVInd[0])] = mesh.tetrahedraNum + 1;
            } else if (tri2Tet.find(Triangle(triVInd[2], triVInd[0], triVInd[1])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[2], triVInd[0], triVInd[1])] = mesh.tetrahedraNum + 1;
            } else if (tri2Tet.find(Triangle(triVInd[2], triVInd[1], triVInd[0])) != tri2Tet.end()) {
                tri2Tet[Triangle(triVInd[2], triVInd[1], triVInd[0])] = mesh.tetrahedraNum + 1;
            } else {
                tri2Tet[Triangle(triVInd[0], triVInd[1], triVInd[2])] = i;
            }
        }
    }

    for (const auto& triI : tri2Tet) {
        const uint64_t& tetId = triI.second;
        const Triangle& triVInd = triI.first;
        if (tetId < mesh.tetrahedraNum) {
            Vector3d vec1 = mesh.vertexes[triVInd[1]] - mesh.vertexes[triVInd[0]];
            Vector3d vec2 = mesh.vertexes[triVInd[2]] - mesh.vertexes[triVInd[0]];
            int id3 = 0;
            for (int i = 0; i < 4; i++) {
                if (mesh.tetrahedras[tetId][i] != triVInd[0]
                    && mesh.tetrahedras[tetId][i] != triVInd[1]
                    && mesh.tetrahedras[tetId][i] != triVInd[2]) {
                    id3 = mesh.tetrahedras[tetId][i];
                    break;
                }
            }

            Vector3d vec3 = mesh.vertexes[id3] - mesh.vertexes[triVInd[0]];
            Vector3d n = vec1.cross(vec2);
            if (n.dot(vec3) < 0) {
                mesh.surface.push_back(Vector4i(triVInd[0], triVInd[1], triVInd[2], tetId));
            } else {
                mesh.surface.push_back(Vector4i(triVInd[0], triVInd[2], triVInd[1], tetId));
            }
        }
    }
    for (const auto& tri : mesh.triangles) {
        mesh.surface.emplace_back(Vector4i(tri[0], tri[1], tri[2], 0));
    }
    vector<bool> flag(mesh.vertexNum, false);
    for (const auto& cTri : mesh.surface) {
        for (int i = 0; i < 3; i++) {
            if (!flag[cTri[i]]) {
                mesh.surfVerts.push_back(cTri[i]);
                flag[cTri[i]] = true;
            }
        }
    }

    set<pair<uint64_t, uint64_t>> SFEdges_set;
    for (const auto& cTri : mesh.surface) {
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
    mesh.surfEdges = vector<pair<uint64_t, uint64_t>>(SFEdges_set.begin(), SFEdges_set.end());
}

} // namespace MeshGeometry
