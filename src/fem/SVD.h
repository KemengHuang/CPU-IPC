#pragma once

#include "Eigen/Eigen"

struct SVDResult2D_double
{
    Eigen::Matrix2d U;
    Eigen::Matrix2d SIGMA;
    Eigen::Matrix2d V;
};

struct SVDResult3D_float
{
    Eigen::Matrix3f U;
    Eigen::Matrix3f SIGMA;
    Eigen::Matrix3f V;
};

struct SVDResult3D_double
{
    Eigen::Matrix3d U;
    Eigen::Matrix3d SIGMA;
    Eigen::Matrix3d V;
};


SVDResult2D_double SingularValueDecomposition2D_double(const Eigen::Matrix2d& F);
SVDResult3D_float SingularValueDecomposition3D_float(const Eigen::Matrix3f& F);
SVDResult3D_double SingularValueDecomposition3D_double(const Eigen::Matrix3d& F);
SVDResult3D_double QRSVD(const Eigen::Matrix3d& F);
