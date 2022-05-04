#pragma once
#include <modus/system/system.hpp>
#include <modus/system/state.hpp>
#include <modus/system/body.hpp>
#include <modus/kinematics/jacobian_floating.hpp>
#include <sophus/se3.hpp>


namespace modus {
namespace examples {

Body* AddSE3Body(System* system);

// se3 twist coordinates, ordering is [x, y, z, w_x, w_y, w_z]
class SE3State : public State {
 protected:
  Body* body_;
  size_t first_dof_index_;
  Eigen::VectorXd q_, qd_, qdd_;
  Sophus::SE3d g_;

 public:
  SE3State(Body* body);

  void SetFirstDofIndex(size_t first) { first_dof_index_ = first; }

  size_t NumDofs() override { return 6; }
  size_t GetFirstDofIndex() override { return first_dof_index_; }

  Eigen::VectorXd GetPositions() override;
  void SetPositions(const Eigen::VectorXd& q) override;

  Eigen::VectorXd GetVelocities() override { return qd_; }
  void SetVelocities(const Eigen::VectorXd& qd) override { qd_ = qd; }

  Eigen::VectorXd GetAccelerations() override { return qdd_; }
  void SetAccelerations(const Eigen::VectorXd& qdd) override { qdd_ = qdd; }
};

}
}