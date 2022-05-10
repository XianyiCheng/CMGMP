#pragma once
#include <vector>
#include <modus/common/eigen.hpp>

namespace modus 
{

using Simplices = std::vector<Eigen::VectorXi>;

using Points = std::vector<Eigen::Vector3d>;

// Convex hull of input points as dxn matrix.
Eigen::MatrixXi ConvexHull(const Eigen::MatrixXd& points, double eps);

// Convex hull with std vector input scheme, will also reorient.
Simplices ConvexHull(const Points& points, double eps);

Simplices ConvertIncidenceMatrixToSimplices(const Eigen::MatrixXi& M);

Simplices ReorderBoundary(const Simplices& simplices, const Points& points);

Simplices ReorientHull(const Simplices& simplices, const Points& points);

}