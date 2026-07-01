#include "collision/BVH.h"
#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

using namespace cipc;

static AABB makePointAABB(const Eigen::Vector3d& p, double eps = 1e-9)
{
    return AABB(p.array() - eps, p.array() + eps);
}

int main()
{
    // Build a BVH over random 3D points and verify that AABB and sphere queries
    // return the same primitive set as a brute-force scan.
    const int n = 1000;
    std::vector<AABB> boxes;
    boxes.reserve(n);
    for (int i = 0; i < n; ++i) {
        Eigen::Vector3d p(
            static_cast<double>(rand()) / RAND_MAX,
            static_cast<double>(rand()) / RAND_MAX,
            static_cast<double>(rand()) / RAND_MAX);
        boxes.push_back(makePointAABB(p));
    }

    BVH bvh;
    bvh.build(boxes);

    const Eigen::Vector3d center(0.5, 0.5, 0.5);
    const double radius = 0.25;

    // Brute force sphere query.
    std::set<int> brute;
    for (int i = 0; i < n; ++i) {
        if ((boxes[i].center() - center).squaredNorm() <= radius * radius) {
            brute.insert(i);
        }
    }

    // BVH query using a conservative box and then exact point-sphere test.
    std::set<int> result;
    AABB queryBox(center.array() - radius, center.array() + radius);
    bvh.queryAABB(queryBox, [&](int pid) {
        if ((boxes[pid].center() - center).squaredNorm() <= radius * radius) {
            result.insert(pid);
        }
    });

    if (result != brute) {
        std::cerr << "BVH sphere query mismatch: got " << result.size()
                  << " expected " << brute.size() << std::endl;
        return 1;
    }

    // Refit test: move all points and rebuild with refit.
    std::vector<AABB> movedBoxes = boxes;
    for (int i = 0; i < n; ++i) {
        movedBoxes[i] = makePointAABB(boxes[i].center() + Eigen::Vector3d(0.1, 0.1, 0.1));
    }
    bvh.refit(movedBoxes);

    const Eigen::Vector3d movedCenter = center + Eigen::Vector3d(0.1, 0.1, 0.1);
    std::set<int> refitResult;
    AABB refitQueryBox(movedCenter.array() - radius, movedCenter.array() + radius);
    bvh.queryAABB(refitQueryBox, [&](int pid) {
        if ((movedBoxes[pid].center() - movedCenter).squaredNorm() <= radius * radius) {
            refitResult.insert(pid);
        }
    });

    std::set<int> refitBrute;
    for (int i = 0; i < n; ++i) {
        if ((movedBoxes[i].center() - movedCenter).squaredNorm() <= radius * radius) {
            refitBrute.insert(i);
        }
    }

    if (refitResult != refitBrute) {
        std::cerr << "BVH refit query mismatch: got " << refitResult.size()
                  << " expected " << refitBrute.size() << std::endl;
        return 1;
    }

    std::cout << "BVH tests passed." << std::endl;
    return 0;
}
