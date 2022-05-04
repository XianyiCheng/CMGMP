#include <modus/geometry/transforms.h>


void modus::triangleFrame(const Eigen::Vector3f& v0,
                         const Eigen::Vector3f& v1,
                         const Eigen::Vector3f& v2,
                         Eigen::Matrix3f& R,
                         Eigen::Vector3f& t) {
    // Translation.
    t = (v0 + v1 + v2)/3.0;
    // Rotation.
    Eigen::Vector3f x = (v1 - v0).normalized();
    Eigen::Vector3f y = (v2 - v0).normalized();
    Eigen::Vector3f z = x.cross(y).normalized();
    y = z.cross(x).normalized();
    R.col(0) = x;
    R.col(1) = y;
    R.col(2) = z;
}

void modus::inverse(Eigen::Matrix3f& R, Eigen::Vector3f& t) {
    R.transposeInPlace();
    t = -R*t;
}