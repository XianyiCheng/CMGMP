#pragma once
#include <modus/kinematics/jacobian.hpp>
#include <modus/system/link.hpp>

namespace modus
{

// Body Jacobian for an unconstrained rigid body.
class FloatingDofsJacobian : public Jacobian {
 public:
  Link* link_;

  FloatingDofsJacobian(Link* link);

  Eigen::RowVectorXd 
  GetGradientAtPointInDirection(const Eigen::Vector3d& point, 
                                const Eigen::Vector3d& direction) override;
};

}