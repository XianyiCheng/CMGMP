#pragma once
#include <modus/modes/graph/contact_mode_graph.hpp>

namespace modus
{

class NodeWrapper : public ContactMode {
 public:
  
};

using NodeWrapperPtr = std::shared_ptr<NodeWrapper>;

class IncidenceGraphWrapper : public ContactModeGraph {
 public:
  
};

using IncidenceGraphWrapperPtr = std::shared_ptr<IncidenceGraphWrapper>;

}