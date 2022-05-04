#pragma once
#include <memory>
#include <modus/common/eigen.hpp>
#include <modus/qpcc/qp.hpp>
#include <qpOASES.hpp>


namespace modus
{

class qpOASESSolver {
 public:
  using real_t = qpOASES::real_t;
  using MatrixXr = Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using VectorXr = Eigen::Matrix<real_t, Eigen::Dynamic, 1, Eigen::ColMajor>;
  using QProblemPtr = std::shared_ptr<qpOASES::QProblem>;

  MatrixXr H_;
  VectorXr g_;
  MatrixXr A_;
  VectorXr lb_;
  VectorXr ub_;
  VectorXr lbA_;
  VectorXr ubA_;
  VectorXr x_;
  QProblemPtr qp_solver_;

  Eigen::VectorXd GetSolution();
  real_t GetCost();
  void Solve(QPPtr qp);
  void InitQP(QPPtr qp);
  void UpdateVectors(QPPtr qp);
};

}