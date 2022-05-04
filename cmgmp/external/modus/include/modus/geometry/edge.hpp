#pragma once
#include <memory>

#include <modus/geometry/halfedge.hpp>

namespace modus
{

class Edge {
public:
    HalfedgePtr halfedge;
    int index;
};

typedef std::shared_ptr<Edge> EdgePtr;

}