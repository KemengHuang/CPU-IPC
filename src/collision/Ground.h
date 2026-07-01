#pragma once

#include "mesh/Mesh.h"
#include <Eigen/Core>

class Ground {
public:
    Eigen::Vector3d normal;
    double D;

    Ground();
    void init(const Eigen::Vector3d& m_normal, const double& d);
    double calculateGapFromObj(const mesh3D& mesh, const int& vId) const;
    void calculateActivateSet(mesh3D& mesh) const;
};
