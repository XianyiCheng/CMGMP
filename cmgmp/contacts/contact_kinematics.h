#ifndef UTILS_H
#define UTILS_H
    #include "../utilities/utilities.h"
#endif

#ifndef SAMPLE_H
#define SAMPLE_H
    #include "../utilities/sample.h"
#endif

#include <math.h>

typedef Eigen::Matrix<double, 7, 1> Vector7d;



Matrix6d contact_jacobian(const Vector3d &position, const Vector3d &normal);

Matrix4d pose2SE3(const Vector7d& x);
Vector7d SE32pose(const Matrix4d& T);
Vector6d compute_rbvel_body(const Vector7d& x, const Vector7d& x_goal);
Vector6d compute_rbvel_spatial(const Vector7d& x, const Vector7d& x_goal);

