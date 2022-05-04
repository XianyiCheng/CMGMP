#pragma once
#include <cstdlib>
#include <modus/common/memory.hpp>
#include <modus/system/aspect.hpp>

namespace modus
{

class Input : public Aspect {//, public std::enable_shared_from_this<Input> {
 public:
  virtual size_t NumInputs() = 0;
  virtual size_t GetFirstInputIndex() = 0;

  virtual Eigen::VectorXd GetInputs() = 0;
  virtual void SetInputs(const Eigen::VectorXd& u) = 0;
};
MODUS_DEFINE_SHARED(Input);

}
