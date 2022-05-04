#pragma once
#include <modus/common/memory.hpp>
#include <modus/common/eigen.hpp>
#include <modus/system/aspect.hpp>


namespace modus
{

class State : public Aspect { //, public std::enable_shared_from_this<State> {
 public:
  virtual size_t NumDofs() = 0;
  virtual size_t GetFirstDofIndex() = 0;

  virtual Eigen::VectorXd GetPositions() = 0;
  virtual void SetPositions(const Eigen::VectorXd& q) = 0;

  virtual Eigen::VectorXd GetVelocities() = 0;
  virtual void SetVelocities(const Eigen::VectorXd& qd) = 0;

  virtual Eigen::VectorXd GetAccelerations() = 0;
  virtual void SetAccelerations(const Eigen::VectorXd& qdd) = 0;
};
MODUS_DEFINE_SHARED(State);

}