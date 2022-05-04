#pragma once
#include <modus/common/memory.hpp>
#include <modus/common/eigen.hpp>
#include <modus/system/aspect.hpp>

namespace modus
{

// Jacobian at a link.
class Jacobian : public Aspect {
 public:
  /**
   * @brief Let d(x) be the position of a point in space along a direction,
   * where x is the body state. This function returns the gradient vector ∂d/∂x.
   *
   * @param point 
   * @param direction 
   * @return Eigen::RowVectorXd
   */
  virtual Eigen::RowVectorXd 
  GetGradientAtPointInDirection(const Eigen::Vector3d& point, 
                                const Eigen::Vector3d& direction) = 0;
};
MODUS_DEFINE_SHARED(Jacobian);

}