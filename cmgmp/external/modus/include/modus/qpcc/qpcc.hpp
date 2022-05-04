#pragma once
#include <memory>
#include <modus/common/eigen.hpp>
#include <modus/qpcc/qp.hpp>
#include <modus/qpcc/cost.hpp>
#include <modus/qpcc/bounds.hpp>
#include <modus/system/system.hpp>


namespace modus
{

// Quadratic program with contact constraints. 
//
// Variable ordering is [q, q̇, τ, λ], where τ are the actuator forces and λ are
// the contact forces.
class QuadraticProgramWithContactConstraints {
 public:
  Costs costs_;
  VariableBounds bounds_;
  System* system_;
  QPPtr qp_;

  void SetSystem(System* system) { system_ = system; }

  void AddCost(CostPtr cost);

  void AddVariableBound(VariableBoundPtr bound);

  // Create quadratic program cost, constraint, and bound matrices.
  QPPtr CreateQuadraticProgram();

  // Update quadratic program's bounds.
  QPPtr UpdateQuadraticProgramBounds(const std::string& cs_mode,
                                     const std::string& ss_mode);
};

using QPCC = QuadraticProgramWithContactConstraints;

using QPCCPtr = std::shared_ptr<QPCC>;

}