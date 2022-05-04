#include <modus/geometry/distance.hpp>
#include <iostream>


float modus::pointTriangleDist2(const Eigen::Vector3d& p,
                               const Eigen::Vector3d& v0,
                               const Eigen::Vector3d& v1,
                               const Eigen::Vector3d& v2) {
    // Calculate coefficients.
    const Eigen::Vector3d e0 = v1 - v0;
    const Eigen::Vector3d e1 = v2 - v0;
    const Eigen::Vector3d vp = v0 - p;
    float a = e0.dot(e0);
    float b = e0.dot(e1);
    float c = e1.dot(e1);
    float d = e0.dot(vp);
    float e = e1.dot(vp);
    float f = vp.dot(vp);

    // Calculate minimum of Q = as^2 + 2bst + ct^2 + 2ds + 2et + f, in scaled
    // coordinates.
    float det = a*c - b*b;
    float s = b*e - c*d;
    float t = b*d - a*e;

    // Handle regions.
    if (s + t <= det) {
        if (s < 0) {
            if (t < 0) {
                // 4
                if (d < 0) {
                    // 5
                    t = 0;
                    s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
                } else {
                    // 3
                    s = 0;
                    t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
                }
            } else {
                // 3
                // F'(t) = 0 when t = -e/c
                s = 0;
                t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
            }
        } else {
            if (t < 0) {
                // 5
                // F'(s) = 0 when s = -d/a
                t = 0;
                s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
            } else {
                // 0
                float invDet = 1.0 / det;
                s *= invDet;
                t *= invDet;
            }
        }
    } else {
        if (s < 0) {
            // 2
            float tmp0 = b + d;
            float tmp1 = c + e;
            if (tmp1 > tmp0) {
                // 1
                float numer = tmp1 - tmp0;
                float denom = a - 2*b + c;
                s = (numer >= denom ? 1 : numer/denom);
                t = 1 - s;
            } else {
                // 3
                s = 0;
                t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : -e/c));
            }
        } else if (t < 0) {
            // 6
            float tmp = a + d;
            if (tmp > 0) {
                // 5
                t = 0;
                s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
            } else {
                // 1
                float numer = c + e - b - d;
                if (numer <= 0) {
                    s = 0;
                } else {
                    float denom = a - 2*b + c;
                    s = (numer >= denom ? 1 : numer/denom);
                }
                t = 1 - s;
            }
        } else {
            // 1
            // F(s) = Q(s, 1 - s) = 
            //      (a - 2b + c)s^2 + (2b - 2c + 2d - 2e)s + c + 2e + f
            // F'(s)/2 = (a - 2b + c)s + b - c + d - e
            // F'(s) = 0 when s = (c + e - b - d)/(a - 2b + c)
            // a - 2b + c = |e0 - e1|^2 > 0
            // so only the sign of b - c + d - e need be considered
            float numer = c + e - b - d;
            if (numer <= 0) {
                s = 0;
            } else {
                float denom = a - 2*b + c;
                s = (numer >= denom ? 1 : numer/denom);
            }
            t = 1 - s;
        }
    }

    return (p - v0 - s*e0 - t*e1).squaredNorm();
}

float modus::pointAABBDist2(const Eigen::Vector3d& p,
                           const Eigen::Vector3d& b_min,
                           const Eigen::Vector3d& b_max) {
    // 
    float cx, cy, cz;
    float ex, ey, ez;
    float vx, vy, vz;
    float rx, ry, rz;

    cx = (b_max[0] + b_min[0])/2.0;
    cy = (b_max[1] + b_min[1])/2.0;
    cz = (b_max[2] + b_min[2])/2.0;

    ex = (b_max[0] - b_min[0])/2.0;
    ey = (b_max[1] - b_min[1])/2.0;
    ez = (b_max[2] - b_min[2])/2.0;

    vx = p[0] - cx;
    vy = p[1] - cy;
    vz = p[2] - cz;

    vx = (abs(vx) <= ex ? vx : (vx > 0 ? ex : -ex));
    vy = (abs(vy) <= ey ? vy : (vy > 0 ? ey : -ey));
    vz = (abs(vz) <= ez ? vz : (vz > 0 ? ez : -ez));

    rx = cx + vx - p[0];
    ry = cy + vy - p[1];
    rz = cz + vz - p[2];

    return rx*rx + ry*ry + rz*rz;
}