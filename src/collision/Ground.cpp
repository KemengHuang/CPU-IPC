#include "collision/Ground.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>

using namespace Eigen;

void Ground::init(const Vector3d& m_normal, const double& d)
{
    normal = m_normal;
    D = d;
}

Ground::Ground()
{
    init(Vector3d(0, 1, 0), -1);
}

double Ground::calculateGapFromObj(const mesh3D& mesh, const int& vId) const
{
    double dist = normal.dot(mesh.vertexes[vId]) - D;
    return dist * dist;
}

void Ground::calculateActivateSet(mesh3D& mesh) const
{
    mesh.Environment_ActiveSet.resize(0);

    tbb::spin_mutex groundMutex;
    tbb::parallel_for(0, (int)mesh.surfVerts.size(), 1, [&](int i) {
        double dis = calculateGapFromObj(mesh, mesh.surfVerts[i]);
        if (dis < mesh.Hhat) {
            tbb::spin_mutex::scoped_lock lock(groundMutex);
            mesh.Environment_ActiveSet.push_back(mesh.surfVerts[i]);
        }
    });
}
