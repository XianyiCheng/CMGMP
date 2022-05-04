#pragma once
#include <modus/common/eigen.hpp>


namespace modus
{

class Dynamics {
 public:
  virtual size_t NumContactForces() = 0;
  virtual size_t NumActuators() = 0;
 
  virtual void GetDynamicEquationsOfMotion(Eigen::MatrixXd& eq_lhs, 
                                           Eigen::VectorXd& eq_rhs) = 0;
  virtual void GetContactForceBounds(Eigen::VectorXd& lb_cf, 
                                     Eigen::VectorXd& ub_cf) = 0;
  virtual void GetActuatorForceBounds(Eigen::VectorXd& lb_cf, 
                                      Eigen::VectorXd& ub_cf) = 0;
};
using DynamicsPtr = std::shared_ptr<Dynamics>;


class DynamicEquationsOfMotionConstraint {
 public:
  double h_;

  void SetTimestep(double h) { h_ = h; };

  virtual size_t NumContactForces() = 0;
  virtual size_t NumActuators() = 0;
  virtual void GetDynamicEquationsOfMotion(Eigen::MatrixXd& eq_lhs, 
                                           Eigen::VectorXd& eq_rhs) = 0;
};
using DynamicEquationsOfMotionConstraintPtr = 
  std::shared_ptr<DynamicEquationsOfMotionConstraint>;
using DynamicEoMConstraintPtr = DynamicEquationsOfMotionConstraintPtr;

}