#pragma once
#include <modus/common/eigen.hpp>
#include <modus/system/aspect.hpp>

namespace modus
{

class GravitoIntertialDynamics : public Aspect
{
 public:
  virtual Eigen::MatrixXd GetMassMatrix() = 0;
  virtual Eigen::MatrixXd GetCoriolisMatrix() = 0;
  virtual Eigen::VectorXd GetPotentialForces() = 0;
};
MODUS_DEFINE_SHARED(GravitoIntertialDynamics);

}