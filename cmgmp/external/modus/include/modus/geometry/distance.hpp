#pragma once
#include <modus/common/eigen.hpp>


namespace modus {

float pointTriangleDist2(const Eigen::Vector3d& p,
                         const Eigen::Vector3d& v0,
                         const Eigen::Vector3d& v1,
                         const Eigen::Vector3d& v2);

float pointAABBDist2(const Eigen::Vector3d& p,
                     const Eigen::Vector3d& b_min,
                     const Eigen::Vector3d& b_max);

}