#ifndef UTILS_H
#define UTILS_H
    #include "utilities.h"
#endif

#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace Eigen;

Vector3d sample_position(const Vector3d &ub, const Vector3d &lb);
Quaterniond generate_unit_quaternion();
Quaterniond sample_rotation(const Vector3d& axis);
Vector3d steer_position(const Vector3d &p0, const Eigen::Vector3d &p1, const double &length);
Quaterniond steer_quaternion(const Quaterniond &q0, const Quaterniond &q1, const double &angle);

void combinations(const std::vector<int> &elems, int req_len, std::vector<std::vector<int>>* elem_combs);
