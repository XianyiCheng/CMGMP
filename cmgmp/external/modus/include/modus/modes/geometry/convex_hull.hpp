#pragma once
#include <modus/common/eigen.hpp>

/**
 * @brief Compute the convex hull of the input points.
 * 
 * @param pts 
 * @param eps 
 * @return Eigen::MatrixXi  vertex-facet incidence matrix
 */
Eigen::MatrixXi convex_hull(const Eigen::MatrixXd& pts, double eps);