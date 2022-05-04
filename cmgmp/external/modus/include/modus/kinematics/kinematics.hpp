#pragma once
#include <modus/common/eigen.hpp>
#include <modus/system/aspect.hpp>
#include <modus/kinematics/jacobian.hpp>

namespace modus
{

// Interface for multiple (possibly articulated) rigid body kinematics.
class Kinematics : public Aspect {
 public:
  virtual Jacobian* GetJacobian(int body_id, int link_id) = 0;
};

}