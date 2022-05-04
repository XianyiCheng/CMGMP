#pragma once
#include <modus/common/eigen.hpp>


Eigen::VectorXd interior_point(const Eigen::MatrixXd& A,
                               const Eigen::VectorXd& b, double eps); //,
                            //    const Eigen::MatrixXd* Aeq=nullptr,
                            //    const Eigen::VectorXd* beq=nullptr);

namespace modus
{
/**
 * @brief Compute an interior point satisfying input cs modes.
 * 
 * @param N         n x d matrix of normal constraints
 * @param cs_mode   contacting-separating mode
 * @param eps       tolerance
 * @return Eigen::VectorXd 
 */
Eigen::VectorXd InteriorPoint(const Eigen::MatrixXd& N, 
                              const std::string& cs_mode,
                              double eps);

/**
 * @brief Compute an interior point satisfying the input cs + ss modes.
 * 
 * @param N         n x d matrix of normal constraints
 * @param T         nk x d matrix of tangent constraints
 * @param cs_mode   contacting-separating mode
 * @param ss_mode   sliding-sticking mode 
 * @param eps       tolerance
 * @return Eigen::VectorXd 
 */
Eigen::VectorXd InteriorPoint(const Eigen::MatrixXd& N,
                              const Eigen::MatrixXd& T,
                              const std::string& cs_mode,
                              const std::string& ss_mode,
                              double eps);
}