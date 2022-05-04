#pragma once
#include <modus/kinematics/jacobian.hpp>
#include <modus/system/link.hpp>

namespace modus
{

class BodyJacobian : public Jacobian {
 public:
  Link* link_;
  Eigen::MatrixXd body_;

  BodyJacobian(Link* link);

  void SetBodyJacobian(const Eigen::MatrixXd& body) { body_ = body; }

  Eigen::RowVectorXd 
  GetGradientAtPointInDirection(const Eigen::Vector3d& point, 
                                const Eigen::Vector3d& direction) override;
};

}