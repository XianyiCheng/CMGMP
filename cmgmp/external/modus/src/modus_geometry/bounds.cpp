#include <modus/geometry/bounds.hpp>
#include <modus/geometry/distance.hpp>


using namespace modus;

Eigen::Vector3d min(const Eigen::Vector3d& p, const Eigen::Vector3d& q) {
    return Eigen::Vector3d(std::min(p[0], q[0]), std::min(p[1], q[1]), 
                           std::min(p[2], q[2]));
}

Eigen::Vector3d max(const Eigen::Vector3d& p, const Eigen::Vector3d& q) {
    return Eigen::Vector3d(std::max(p[0], q[0]), std::max(p[1], q[1]), 
                           std::max(p[2], q[2]));
}

Bounds3 modus::unionBB(const Bounds3& b0, const Bounds3& b1) {
    Bounds3 ret;
    ret.pMin = min(b0.pMin, b1.pMin);
    ret.pMax = max(b0.pMax, b1.pMax);
    return ret;
}

Bounds3 modus::unionBP(const Bounds3& b, const Eigen::Vector3d& p) {
    Bounds3 ret;
    ret.pMin = min(b.pMin, p);
    ret.pMax = max(b.pMax, p);
    return ret;
}

Bounds3::Bounds3() {
    float minNum = -std::numeric_limits<float>::infinity();
    float maxNum =  std::numeric_limits<float>::infinity();
    pMin = Eigen::Vector3d(maxNum, maxNum, maxNum);
    pMax = Eigen::Vector3d(minNum, minNum, minNum);
}

Bounds3::Bounds3(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) 
    : pMin(std::min(p1[0], p2[0]), std::min(p1[1], p2[1]), 
           std::min(p1[2], p2[2])),
      pMax(std::max(p1[0], p2[0]), std::max(p1[1], p2[1]), 
           std::max(p1[2], p2[2])) {
}

Eigen::Vector3d Bounds3::diagonal() {
    return pMax - pMin;
}

Eigen::Vector3d Bounds3::center() {
    return (pMax + pMin)/2.0;
}

float Bounds3::surfaceArea() {
    const Eigen::Vector3d d = diagonal();
    return 2.0 * (d[0]*d[1] + d[0]*d[2] + d[1]*d[2]);
}

float Bounds3::volume() {
    const Eigen::Vector3d d = diagonal();
    return d[0]*d[1]*d[2];
}

int Bounds3::maxExtent() {
    const Eigen::Vector3d d = diagonal();
    if ((d[0] > d[1]) && (d[0] > d[2])) {
        return 0;
    } else if (d[1] > d[2]) {
        return 1;
    } else {
        return 2;
    }
}

Eigen::Vector3d Bounds3::offset(const Eigen::Vector3d& p) const {
    Eigen::Vector3d o = p - pMin;
    if (pMax[0] > pMin[0]) {
        o[0] /= pMax[0]-pMin[0];
    }
    if (pMax[1] > pMin[1]) {
        o[1] /= pMax[1]-pMin[1];
    }
    if (pMax[2] > pMin[2]) {
        o[2] /= pMax[2]-pMin[2];
    }
    return o;
}

bool Bounds3::intersectRay(RayPtr ray) {

}

float Bounds3::distance2(const Eigen::Vector3d& p) {
    return pointAABBDist2(p, pMin, pMax);
}

std::ostream& operator<<(std::ostream& os, const Bounds3& b) {
    os << "bounds:\n" << b.pMin.transpose() << "\n" << b.pMax.transpose();
    return os;
}