#pragma once
#include <modus/common/eigen.hpp>
#include <modus/system/system.hpp>
#include <modus/system/aspect.hpp>


namespace modus
{

class ContinuousDynamics : public Aspect {
 protected:
  System* system_;

 public:
  ContinuousDynamics(System* system);

  Eigen::MatrixXd GetMassMatrix();
  Eigen::MatrixXd GetCoriolisMatrix();
  Eigen::VectorXd GetPotentialForces();
  Eigen::MatrixXd GetActuatorMatrix();
};
MODUS_DEFINE_SHARED(ContinuousDynamics);

}