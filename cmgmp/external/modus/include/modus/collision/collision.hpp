#pragma once
#include <modus/common/memory.hpp>
#include <modus/collision/contact.hpp>
#include <modus/system/aspect.hpp>

namespace modus
{

class Collision : public Aspect {
 public:
  virtual Contacts GetContactPoints() = 0;
};
MODUS_DEFINE_SHARED(Collision);

}