#pragma once
#include <vector>
#include <map>
#include <modus/common/memory.hpp>
#include <modus/common/eigen.hpp>
#include <modus/qpcc/qpcc.hpp>


namespace modus
{

class QPCCSolver {
 public:
  virtual void Solve(QPCCPtr qpcc) = 0;

  virtual Eigen::VectorXd GetSolution() = 0;
};

}