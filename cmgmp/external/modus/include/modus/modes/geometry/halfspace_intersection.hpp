#pragma once
#include <modus/modes/geometry/incidence_graph.hpp>


/**
 * @brief Compute the incidence graph of the halfspace intersection A - b<= 0.
 *
 * @param A 
 * @param b 
 * @param eps 
 * @param int_pt
 * @return IncidenceGraph* 
 */
IncidenceGraph* halfspace_intersection(const Eigen::MatrixXd& A,
                                       const Eigen::VectorXd& b,
                                       double eps,
                                       const Eigen::VectorXd* int_pt=nullptr);