#pragma once
#include <memory>
#include <vector>


namespace modus
{

class ContactMode;
using ContactModePtr = std::shared_ptr<ContactMode>;

class ContactMode {
 public:
  std::string GetCSMode();
  std::string GetSSMode();
  std::vector<ContactModePtr> AdjacentModes();
};

class ContactModeGraph {
 public:
  std::vector<ContactModePtr> GetVertices();
  ContactModePtr GetRoot();
};

}