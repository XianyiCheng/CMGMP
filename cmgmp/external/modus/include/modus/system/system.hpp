#pragma once
#include <memory>
#include <vector>
#include <modus/common/memory.hpp>
#include <modus/common/eigen.hpp>
#include <modus/system/composite.hpp>
#include <modus/system/body.hpp>
#include <modus/system/link.hpp>
#include <modus/system/state.hpp>
#include <modus/system/input.hpp>


namespace modus
{

class System : public Composite {
 protected:
  std::vector<BodyPtr> bodies_;

 public:
  System();

  size_t NumBodies();

  Body* GetBody(int body_id);

  Body* CreateBody();
};
MODUS_DEFINE_SHARED(System);

class SystemState : public State {
 protected:
  System* system_;

 public:
  SystemState(System* system);

  size_t NumDofs() override;
  size_t GetFirstDofIndex() override;

  Eigen::VectorXd GetPositions() override;
  void SetPositions(const Eigen::VectorXd& q) override;

  Eigen::VectorXd GetVelocities() override;
  void SetVelocities(const Eigen::VectorXd& qd) override;

  Eigen::VectorXd GetAccelerations() override;
  void SetAccelerations(const Eigen::VectorXd& qdd) override;
};
MODUS_DEFINE_SHARED(SystemState);

class SystemInput : public Input {
 protected:
  System* system_;

 public:
  SystemInput(System* system);

  size_t NumInputs() override;
  size_t GetFirstInputIndex() override;

  Eigen::VectorXd GetInputs() override;
  void SetInputs(const Eigen::VectorXd& u) override;
};
MODUS_DEFINE_SHARED(SystemInput);

}