#pragma once
#include <memory>
#include <vector>
#include <modus/common/memory.hpp>
#include <modus/common/eigen.hpp>
#include <modus/system/composite.hpp>
#include <modus/system/link.hpp>


namespace modus
{

class Body : public Composite {
 protected:
  int index_;

 public:
  std::vector<LinkPtr> links_;

  Body();

  void SetIndex(int index) { index_ = index; }
  int GetIndex();

  size_t NumLinks();

  Link* GetRootLink();
  Link* GetLink(int link_id);

  Link* CreateLink();
};
MODUS_DEFINE_SHARED(Body);

}