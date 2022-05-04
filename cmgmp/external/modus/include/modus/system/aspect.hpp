#pragma once
#include <modus/common/memory.hpp>


namespace modus
{

class Aspect : public std::enable_shared_from_this<Aspect> {
 public:
  virtual ~Aspect() = default;
};
MODUS_DEFINE_SHARED(Aspect);

}