#pragma once
#include <memory>
#include <modus/common/eigen.hpp>


namespace modus
{

// Quadratic program.
struct QuadraticProgram {
  // Quadratic cost: 
  //  1/2xᵀHx + xᵀg
  Eigen::MatrixXd H;
  Eigen::VectorXd g;
  
  // Constraint matrix and constraint bounds:
  //  lbA ≤ Ax ≤ ubA
  Eigen::MatrixXd A;
  Eigen::VectorXd lbA;
  Eigen::VectorXd ubA;

  // Variable bounds:
  //  lb ≤ x ≤ ub
  Eigen::VectorXd lb;
  Eigen::VectorXd ub;

  QuadraticProgram(int n_v, int n_c);
};

using QP = QuadraticProgram;
using QPPtr = std::shared_ptr<QP>;

std::ostream& operator<<(std::ostream& os, const QP& qp);

}