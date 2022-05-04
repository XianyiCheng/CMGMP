#pragma once

#include <Eigen/Dense>
#include <math.h>
#include <iostream>

#define ZERO_SOLUTION 0
#define ONE_SOLUTION 1
#define TWO_SOLUTIONS 2
#define INF_SOLUTION 3

using namespace std;
class TwoLinkKinematics {
  public:
    TwoLinkKinematics(double l1, double l2);
    Eigen::Vector2d fk(Eigen::Vector2d joints);
    std::tuple<Eigen::Vector2d, Eigen::Vector2d, int> ik(Eigen::Vector2d pose);
    Eigen::Matrix2d jacobian(Eigen::Vector2d joints);
    bool inJointLimits(Eigen::Vector2d joints);
  private:
    double l1_, l2_;
};
