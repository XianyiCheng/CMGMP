#pragma once
#include <modus/common/memory.hpp>
#include <modus/common/eigen.hpp>
#include <modus/system/aspect.hpp>

namespace Eigen {
  using Matrix6d = Matrix<double, 6, 6>;
}

namespace modus
{

class MassMatrix : public Aspect
{
 protected:
  double mass_;
  Eigen::Matrix3d inertia_;

 public:
  double GetMass() { return mass_; }

  Eigen::Matrix3d GetInertia() { return inertia_; }

  // TODO Move to CPP
  Eigen::Matrix6d GetMassMatrix() {
    Eigen::Matrix6d M;
    M.setIdentity();
    M.topRows(3) *= mass_;
    M.bottomRightCorner(3, 3) = inertia_;
    return M;
  }
};
MODUS_DEFINE_SHARED(MassMatrix);

}