#pragma once

#include <string>
#include "Eigen/Eigen"

class mesh3D;

namespace MeshIO {

bool loadTetrahedraMesh(mesh3D& mesh, const std::string& filename, double scale, Eigen::Vector3d offset);
bool loadTriangleMesh(mesh3D& mesh, const std::string& filename, double scale, Eigen::Vector3d offset, int type = 0);
bool outputTetrahedraMesh(const mesh3D& mesh, const std::string& filename);
bool outputTetTempData(const mesh3D& mesh);
bool loadTetTempData(mesh3D& mesh);

} // namespace MeshIO
