#pragma once
#include <modus/common/eigen.hpp>
#include <modus/system/aspect.hpp>

namespace modus
{

class ActuatorDynamics : public Aspect
{
 public:
  virtual Eigen::MatrixXd GetActuatorMatrix() = 0;
};
MODUS_DEFINE_SHARED(ActuatorDynamics);

}