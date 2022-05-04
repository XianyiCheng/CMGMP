#pragma once
#include <modus/system/system.hpp>
#include <modus/integration/integration.hpp>


namespace modus {

class ForwardEuler : public Integration {
 protected:
  double  timestep_;
  System* system_;

 public:
  ForwardEuler(System* system, double timestep=0.01);

  void SetTimestep(double h) { timestep_ = h; }

  void GetODE(Eigen::MatrixXd& lhs, Eigen::VectorXd& rhs) override;
};

}