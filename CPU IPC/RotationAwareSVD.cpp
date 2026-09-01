#include "RotationAwareSVD.h"

#include "ImplicitQR3x3SVD.h"

RotationAwareSVDResult computeRotationAwareSVD(const Eigen::Matrix3d& matrix)
{
    RotationAwareSVDResult result;
    Eigen::Vector3d singularValues;
    JIXIE::singularValueDecomposition(
        matrix,
        result.leftSingularVectors,
        singularValues,
        result.rightSingularVectors);
    result.singularValueMatrix = singularValues.asDiagonal();
    return result;
}
