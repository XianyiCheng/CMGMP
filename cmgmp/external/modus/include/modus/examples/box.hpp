#pragma once
#include <modus/system/system.hpp>
#include <modus/system/body.hpp>
#include <modus/system/state.hpp>
#include <modus/dynamics/gravito_intertial_dynamics.hpp>
#include <modus/dynamics/actuator_dynamics.hpp>
#include <modus/dynamics/mass_matrix.hpp>


namespace modus {
namespace examples {

// Add basic 1-1-1 box to the system. This includes the box state + dynamics +
// actuation. 
Body* AddBox(System* system, bool actuated=false);

class BoxMassMatrix : public MassMatrix {
 public:
  BoxMassMatrix(double mass, double x=1.0, double y=1.0, double z=1.0);
};

// Quasi-static dynamics of a box with respect to the body frame.
class BoxQuasiStatics : public GravitoIntertialDynamics {
 protected:
  Body* box_;

 public:
  BoxQuasiStatics(Body* box);

  Eigen::MatrixXd GetMassMatrix() override;
  Eigen::MatrixXd GetCoriolisMatrix() override;
  Eigen::VectorXd GetPotentialForces() override;
};

// 6-Dof controls for an actuated box.
class BoxControls : public Input {
 public:
  Eigen::VectorXd u_;
  size_t first_input_index_;

  void SetFirstInputIndex(size_t first);

  size_t NumInputs() override;
  size_t GetFirstInputIndex() override;

  Eigen::VectorXd GetInputs() override;
  void SetInputs(const Eigen::VectorXd& u) override;
};

// Actuated dynamics of a box with respect to the body frame.
class BoxActuator: public ActuatorDynamics {
 public:
  BoxActuator(Body* box) {}

  Eigen::MatrixXd GetActuatorMatrix() override;
};

}
}