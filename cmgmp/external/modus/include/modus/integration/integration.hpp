#pragma once
#include <modus/common/eigen.hpp>
#include <modus/system/aspect.hpp>


namespace modus
{

class Integration : public Aspect {
 public:
  virtual void GetODE(Eigen::MatrixXd& lhs, Eigen::VectorXd& rhs) = 0;
};
MODUS_DEFINE_SHARED(Integration);

}