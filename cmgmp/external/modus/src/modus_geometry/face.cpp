#include <iostream>
#include <modus/geometry/face.hpp>
#include <modus/geometry/vertex.hpp>
#include <modus/geometry/distance.hpp>


using namespace modus;

Face::Face(bool isBoundary) : boundary(isBoundary) {
    this->halfedge = nullptr;
}

Bounds3 Face::worldBound() {
    const Eigen::Vector3d& p0 = halfedge->vertex->position;
    const Eigen::Vector3d& p1 = halfedge->next->vertex->position;
    const Eigen::Vector3d& p2 = halfedge->next->next->vertex->position;
    return unionBP(Bounds3(p0, p1), p2);
}

float Face::distance2(const Eigen::Vector3d& p) {
    const Eigen::Vector3d& p0 = halfedge->vertex->position;
    const Eigen::Vector3d& p1 = halfedge->next->vertex->position;
    const Eigen::Vector3d& p2 = halfedge->next->next->vertex->position;
    return pointTriangleDist2(p, p0, p1, p2);
}

int Face::degree() {
    return 0;
}

bool Face::isBoundary() {
    return boundary;
}

Eigen::Vector3d Face::position() {
    Eigen::Vector3d pos;
    position(pos);
    return pos;
}

void Face::position(Eigen::Vector3d& pos) {
    pos.setZero();
    HalfedgePtr h = halfedge;
    int n = 0;
    while (true) {
        pos += h->vertex->position;
        h = h->next;
        n += 1;
        if (h == halfedge) {
            break;
        }
    }
    assert(n == 3);
    pos /= n;
}

Eigen::Vector3d Face::normal() {
    Eigen::Vector3d n;
    n.setZero();
    HalfedgePtr h = halfedge;
    while (true) {
        const Eigen::Vector3d& pi = h->vertex->position;
        const Eigen::Vector3d& pj = h->next->vertex->position;
        n += pi.cross(pj);
        h = h->next;
        if (h == halfedge) {
            break;
        }
    }
    n.normalize();
    return n;
}