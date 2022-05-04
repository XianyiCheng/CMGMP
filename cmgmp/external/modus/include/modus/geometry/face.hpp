#pragma once
#include <memory>
#include <modus/common/eigen.hpp>
// #include <modus/geometry/shape.hpp>
#include <modus/geometry/bounds.hpp>
#include <modus/geometry/halfedge.hpp>


namespace modus
{

class Face {
public:
    HalfedgePtr halfedge;
    bool boundary;
    int index;

    Face(bool isBoundary = false);

    Bounds3 worldBound();
    float distance2(const Eigen::Vector3d& point);

    int degree();
    double area();
    Eigen::Vector3d normal();
    bool isBoundary();
    Eigen::Vector3d position();
    void position(Eigen::Vector3d& pos);
};

typedef std::shared_ptr<Face> FacePtr;

}