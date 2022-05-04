#pragma once
#include <modus/kinematics/jacobian.hpp>
#include <modus/system/link.hpp>

namespace modus
{

/**
 * @brief A geometric jacobian. The linear jacobian is in world frame. The
 * angular jacobian 
 *
 */
class GeometricJacobian : public Jacobian {
 public:
  Link* link_;

  Eigen::MatrixXd linear_;
  Eigen::MatrixXd rotational_;

  GeometricJacobian(Link* link);

  void SetLinearJacobian(const Eigen::MatrixXd& linear);
  void SetRotationalJacobian(const Eigen::MatrixXd& rotational);

  Eigen::RowVectorXd 
  GetGradientAtPointInDirection(const Eigen::Vector3d& point, 
                                const Eigen::Vector3d& direction) override;
};

}