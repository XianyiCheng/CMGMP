#pragma once
#include <memory>
#include <map>
#include <modus/common/memory.hpp>
#include <modus/common/eigen.hpp>
#include <modus/system/composite.hpp>


namespace modus
{

class Link : public Composite {
 protected:
  int index_;

 public:
  void SetIndex(int index) { index_ = index; }
  int GetIndex() { return index_; }
};
MODUS_DEFINE_SHARED(Link);

}