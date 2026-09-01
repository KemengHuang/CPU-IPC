#pragma once

#include <Eigen/Core>

struct RotationAwareSVDResult {
    Eigen::Matrix3d leftSingularVectors;
    Eigen::Matrix3d singularValueMatrix;
    Eigen::Matrix3d rightSingularVectors;
};

// Computes A = U * Sigma * V^T while keeping U and V proper rotations.
// Singular values are ordered by decreasing magnitude; only the final value
// may be negative, which preserves the sign of det(A) for inverted elements.
RotationAwareSVDResult computeRotationAwareSVD(const Eigen::Matrix3d& matrix);
