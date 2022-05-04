#pragma once
#include <memory>
#include <modus/common/eigen.hpp>


namespace modus
{

class VariableBound {
 public:
  virtual void Apply(Eigen::VectorXd& lb, Eigen::VectorXd& ub) = 0;
};

using VariableBoundPtr = std::shared_ptr<VariableBound>;
using VariableBounds = std::vector<VariableBoundPtr>;

}