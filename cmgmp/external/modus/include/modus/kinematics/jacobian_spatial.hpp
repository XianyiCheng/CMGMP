#pragma once
#include <modus/kinematics/jacobian.hpp>
#include <modus/system/link.hpp>

namespace modus
{

class SpatialJacobian : public Jacobian {
 public:
  Eigen::MatrixXd spatial_;

  void SetSpatialJacobian(const Eigen::MatrixXd& spatial) { spatial_ = spatial; }

  Eigen::RowVectorXd 
  GetGradientAtPointInDirection(const Eigen::Vector3d& point, 
                                const Eigen::Vector3d& direction) override;
};

}