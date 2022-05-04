#pragma once
#include <modus/modes/geometry/incidence_graph.hpp>


std::vector<int> 
initial_hyperplanes(Eigen::MatrixXd& A, Eigen::VectorXd& b, double eps);

/**
 * @brief Initialize a hyperplane arrangement.
 * 
 * @param A 
 * @param b 
 * @param eps 
 * @return IncidenceGraph* 
 */
IncidenceGraph* initial_arrangement(const Eigen::MatrixXd& A, 
                                    const Eigen::VectorXd& b, 
                                    double eps);

/**
 * @brief Preprocess hyperplanes for partial arrangement. This function does the
 * following:
 *      - Remove duplicate hyperplanes.
 *      - Reorder hyperplanes to be linearly independent in top d rows.
 *
 * @param A 
 * @param b 
 * @param sides 
 * @param eps 
 */
void partial_preprocess_hyperplanes(Eigen::MatrixXd& A, 
                                    Eigen::VectorXd& b,
                                    Eigen::VectorXi& sides, double eps);

/**
 * @brief Initialize a partial hyperplane arrangement with enumeration.
 *
 * By convention, when sides[i] = 2, then we add the {0,-} sides of hyperplane
 * i, and when sides[i] = 3, we add the {+,0,-} sides. Other values are not
 * supported.
 *
 * @param A     dxd matrix
 * @param b     dx1 vector
 * @param sides dx1 vector
 * @param eps   tolerance
 * @return IncidenceGraph*  
 */
IncidenceGraph* partial_initial_arrangement(const Eigen::MatrixXd& A,
                                            const Eigen::VectorXd& b,
                                            const Eigen::VectorXi& sides,
                                            double eps);

/**
 * @brief Initialize a partial hyperplane arrangement with convex hull.
 *
 * @param A     dxd matrix
 * @param b     dx1 vector
 * @param sides dx1 vector
 * @param eps   tolerance
 * @return IncidenceGraph*  
 */
IncidenceGraph* partial_initial_convex(const Eigen::MatrixXd& A,
                                       const Eigen::VectorXd& b,
                                       double eps);

void increment_arrangement(Eigen::VectorXd a, double b, 
                           IncidenceGraph* I, double eps);