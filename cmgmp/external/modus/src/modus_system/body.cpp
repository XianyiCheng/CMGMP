#include <modus/system/body.hpp>

using namespace modus;


modus::Body::Body() {
  // Add root link.
  CreateLink();
}

int modus::Body::GetIndex() { 
  return index_;
}

size_t modus::Body::NumLinks() {
  return links_.size();
}

Link* modus::Body::GetRootLink() {
  return links_[0].get();
}

Link* modus::Body::GetLink(int link_id) {
  return links_[link_id].get();
}

Link* modus::Body::CreateLink() {
  std::unique_ptr<Link> link = std::make_unique<Link>();
  link->SetIndex(links_.size());
  Link* ret = link.get();
  links_.push_back(std::move(link));
  return ret;
}