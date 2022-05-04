#pragma once
#include <memory>
#include <vector>
#include <modus/common/eigen.hpp>
#include <modus/system/system.hpp>


namespace modus
{

// Quadratic cost: 
//  1/2xᵀHx + xᵀg
class Cost {
 public:
  virtual void AddToObjective(Eigen::MatrixXd& H, Eigen::VectorXd& g) = 0;
};

using CostPtr = std::shared_ptr<Cost>;
using Costs = std::vector<CostPtr>;

class VelocityCost : public Cost {
 protected:
  System*         system_;
  State*          state_;
  Eigen::VectorXd target_velocity_;
  double          scale_; // 

 public:
  VelocityCost(State* state, System* system, const Eigen::VectorXd& target);

  void SetScale(double scale) { scale_ = scale; }

  void AddToObjective(Eigen::MatrixXd& H, Eigen::VectorXd& g);
};

class ContactForceCost : public Cost {
 protected:
  System* system_;
  double  scale_;

 public:
  ContactForceCost(System* system);

  void SetScale(double scale) { scale_ = scale; }

  void AddToObjective(Eigen::MatrixXd& H, Eigen::VectorXd& g);
};

class InputCost : public Cost {
 protected:
  System* system_;
  Input*  input_;
  double  scale_;

 public:
  InputCost(Input* input, System* system);

  void SetScale(double scale) { scale_ = scale; }

  void AddToObjective(Eigen::MatrixXd& H, Eigen::VectorXd& g);
};

}