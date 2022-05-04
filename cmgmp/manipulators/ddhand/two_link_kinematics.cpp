#include "two_link_kinematics.h"

double wrap_to_pi(double value)
{
  if (value > M_PI)
    value -= 2 * M_PI;
  if (value < -M_PI)
    value += 2 * M_PI;
  return value;
}

double deg2rad (double degrees) {
    return degrees * M_PI/180.;
}
double rad2deg (double radians) {
    return radians * 180./M_PI;
}

TwoLinkKinematics::TwoLinkKinematics(double l1, double l2) : l1_(l1), l2_(l2)
{
}

Eigen::Vector2d TwoLinkKinematics::fk(Eigen::Vector2d joints)
{
  double x1 = l1_ * cos(joints[0]);
  double y1 = l1_ * sin(joints[0]);
  double x2 = x1 + l2_ * cos(joints[1]);
  double y2 = y1 + l2_ * sin(joints[1]);

  Eigen::Vector2d pose(x2, y2);
  return pose;
}

std::tuple<Eigen::Vector2d, Eigen::Vector2d, int> TwoLinkKinematics::ik(Eigen::Vector2d pose)
{

  double l1l2sq = pow(l1_, 2) + pow(l2_, 2);
  double c_2 = (pose.squaredNorm() - l1l2sq) / (2 * l1_ * l2_);
  Eigen::Vector2d joints;
  if (abs(c_2) > 1)
  { //No solution
    joints << NAN, NAN;
    return std::make_tuple(joints, joints, ZERO_SOLUTION);
  }
  if (c_2 == 1)
  { // One solution
    joints << atan2(pose[1], pose[0]), 0;
    return std::make_tuple(joints, joints, ONE_SOLUTION);
  }
  if (c_2 == -1)
  {
    if (pose.norm() > 0)
    { //One Solution
      joints << atan2(pose[1], pose[0]), M_PI;
      return std::make_tuple(joints, joints, ONE_SOLUTION);
    }
    else
    { //Inf solutions
      joints << 0, M_PI;
      return std::make_tuple(joints, joints, INF_SOLUTION);
    }
  }
  double x = pose[0];
  double y = pose[1];
  double sigma1 = sqrt(2 * l1l2sq * pose.squaredNorm() - pow(pow(l1_, 2) - pow(l2_, 2), 2) - pow(pose.squaredNorm(), 2));
  double q11 = 2 * atan2(2 * l1_ * y + sigma1, pow(l1_, 2) + 2 * l1_ * x - pow(l2_, 2) + pose.squaredNorm());
  double q12 = 2 * atan2(2 * l1_ * y - sigma1, pow(l1_, 2) + 2 * l1_ * x - pow(l2_, 2) + pose.squaredNorm());
  double q2 = 2 * atan2(sqrt((pose.squaredNorm() - pow(l1_ - l2_, 2)) * (pow(l1_+l2_,2) - pose.squaredNorm())), pose.squaredNorm() - pow(l1_ - l2_, 2));
  double q21 = - q2;
  double q22 =   q2;
  q11 = wrap_to_pi(q11);
  q12 = wrap_to_pi(q12);
  q21 = wrap_to_pi(q21);
  q22 = wrap_to_pi(q22);
  Eigen::Vector2d joints1(q11, q11 + q21);
  Eigen::Vector2d joints2(q12, q12 + q22);
  return std::make_tuple(joints1, joints2, TWO_SOLUTIONS);
}

Eigen::Matrix2d TwoLinkKinematics::jacobian(Eigen::Vector2d joints)
{
  Eigen::Matrix2d jacobian;
  jacobian << -l1_ * sin(joints[0]),
              -l2_ * sin(joints[1]),
              l1_ * cos(joints[0]),
              l2_ * cos(joints[1]);
  return jacobian;
}
bool TwoLinkKinematics::inJointLimits(Eigen::Vector2d joints){
  const double w1 = 8e-3;
  const double w2 = 9e-3;
  const double w3 = 8e-3;
  const double w4 = 8e-3;

  const double d1 = l1_*sin(joints[0]-joints[1]);
  const double d2 = 15e-3*sin(joints[1] + deg2rad(60));

  if (abs(d1)< (w1+w2)/2)
    return false;
  if (d2 < (w3+w4)/2)
    return false;
  if (joints[0] < -deg2rad(140) && joints[1] > -deg2rad(70))
    return false;

  return true;




}
