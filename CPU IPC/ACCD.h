#pragma once
#include "Eigen/Eigen"

double point_triangle_ccd(
    const Eigen::Vector3d& _p,
    const Eigen::Vector3d& _t0,
    const Eigen::Vector3d& _t1,
    const Eigen::Vector3d& _t2,
    const Eigen::Vector3d& _dp,
    const Eigen::Vector3d& _dt0,
    const Eigen::Vector3d& _dt1,
    const Eigen::Vector3d& _dt2,
    double eta, double thickness);

double edge_edge_ccd(
    const Eigen::Vector3d& _ea0,
    const Eigen::Vector3d& _ea1,
    const Eigen::Vector3d& _eb0,
    const Eigen::Vector3d& _eb1,
    const Eigen::Vector3d& _dea0,
    const Eigen::Vector3d& _dea1,
    const Eigen::Vector3d& _deb0,
    const Eigen::Vector3d& _deb1,
    double eta, double thickness);

bool point_triangle_ccd_broadphase(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& t0,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2,
    const Eigen::Vector3d& dp,
    const Eigen::Vector3d& dt0,
    const Eigen::Vector3d& dt1,
    const Eigen::Vector3d& dt2,
    double dist);

bool edge_edge_ccd_broadphase(
    const Eigen::Vector3d& ea0,
    const Eigen::Vector3d& ea1,
    const Eigen::Vector3d& eb0,
    const Eigen::Vector3d& eb1,
    const Eigen::Vector3d& dea0,
    const Eigen::Vector3d& dea1,
    const Eigen::Vector3d& deb0,
    const Eigen::Vector3d& deb1,
    double dist);
