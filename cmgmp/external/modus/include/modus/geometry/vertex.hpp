#pragma once
#include <memory>
#include <modus/common/eigen.hpp>

#include <modus/geometry/halfedge.hpp>

namespace modus
{

class Vertex {
public:
    HalfedgePtr     halfedge;
    Eigen::Vector3d position;
    int index;

    Vertex();

    Eigen::Vector3d normal();
};

typedef std::shared_ptr<Vertex> VertexPtr;

}