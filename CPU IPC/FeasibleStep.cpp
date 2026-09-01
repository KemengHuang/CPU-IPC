#include "FeasibleStep.h"

#include <tbb/parallel_for.h>

#include <algorithm>
#include <stdexcept>

void limitStepByGround(
    const mesh3D& mesh,
    const Ground& ground,
    const std::vector<Eigen::Vector3d>& searchDirection,
    double slackness,
    double& stepSize)
{
    if (mesh.surfVerts.empty()) {
        return;
    }
    if (searchDirection.size() != mesh.vertexes.size()
        || slackness <= 0.0 || slackness > 1.0) {
        throw std::invalid_argument("invalid ground feasible-step input");
    }

    Eigen::VectorXd maximumSteps(mesh.surfVerts.size());
    tbb::parallel_for(0, static_cast<int>(mesh.surfVerts.size()), 1,
        [&](int surfaceVertex) {
            maximumSteps[surfaceVertex] = 1.0;
            const int vertex = static_cast<int>(mesh.surfVerts[surfaceVertex]);
            const double normalMotion = ground.normal.dot(-searchDirection[vertex]);
            if (normalMotion >= 0.0) {
                return;
            }

            const double distance = ground.normal.dot(mesh.vertexes[vertex]) - ground.D;
            const double candidate = -distance / normalMotion * slackness;
            if (candidate >= 0.0) {
                maximumSteps[surfaceVertex] = candidate;
            }
        });
    stepSize = (std::min)(stepSize, maximumSteps.minCoeff());
}
