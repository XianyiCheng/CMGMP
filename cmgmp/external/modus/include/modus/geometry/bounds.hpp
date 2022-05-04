#pragma once
#include <iostream>
#include <memory>
#include <algorithm>
#include <modus/common/eigen.hpp>
#include <modus/geometry/ray.hpp>


namespace modus {

class Bounds3 {
public:
    Eigen::Vector3d pMin;
    Eigen::Vector3d pMax;

    Bounds3();
    Bounds3(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);

    Eigen::Vector3d diagonal();
    Eigen::Vector3d center();
    float surfaceArea();
    float volume();
    int maxExtent();
    Eigen::Vector3d offset(const Eigen::Vector3d& p) const;

    bool intersectRay(RayPtr ray);
    float distance2(const Eigen::Vector3d& p);

    friend std::ostream& operator<<(std::ostream& os, const Bounds3& b);
};

Bounds3 unionBB(const Bounds3& b0, const Bounds3& b1);

Bounds3 unionBP(const Bounds3& b, const Eigen::Vector3d& p);

typedef std::shared_ptr<Bounds3> Bounds3Ptr;

}

