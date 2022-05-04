#include <modus/geometry/halfedge.hpp>
#include <modus/geometry/vertex.hpp>
#include <modus/geometry/face.hpp>


using namespace modus;

Halfedge::Halfedge() : index(0) {
    this->next = nullptr;
    this->twin = nullptr;
    this->vertex = nullptr;
    this->edge = nullptr;
    this->face = nullptr;
}

bool Halfedge::isBoundary() {
    return face->isBoundary();
}


void Halfedge::setNeighbors(HalfedgePtr n, HalfedgePtr t, VertexPtr v, 
                            EdgePtr e, FacePtr f) {
    this->next = n;
    this->twin = t;
    this->vertex = v;
    this->edge = e;
    this->face = f;
}

Eigen::Vector3d Halfedge::vector() {
    return this->next->vertex->position - this->vertex->position;
}