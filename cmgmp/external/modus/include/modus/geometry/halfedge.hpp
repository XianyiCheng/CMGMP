#pragma once
#include <memory>
#include <modus/common/eigen.hpp>

namespace modus
{

class Halfedge;
class Vertex;
class Edge;
class Face;

typedef std::shared_ptr<Halfedge> HalfedgePtr;
typedef std::shared_ptr<Vertex> VertexPtr;
typedef std::shared_ptr<Edge> EdgePtr;
typedef std::shared_ptr<Face> FacePtr;

class Halfedge {
public:
    HalfedgePtr next;
    HalfedgePtr twin;
    VertexPtr   vertex;
    EdgePtr     edge;
    FacePtr     face;
    int index;

    Halfedge();

    bool isBoundary();

    Eigen::Vector3d vector();

    void setNeighbors(HalfedgePtr next, HalfedgePtr twin, VertexPtr vertex, 
                      EdgePtr edge, FacePtr face);
};

typedef std::shared_ptr<Halfedge> HalfedgePtr;

}