#include <modus/kinematics/jacobian_geometric.hpp>
#include <modus/system/transform.hpp>
#include <modus/common/linear_algebra.hpp>

using namespace modus;


modus::GeometricJacobian::GeometricJacobian(Link* link) {
  link_ = link;
}

void modus::GeometricJacobian::SetLinearJacobian(const Eigen::MatrixXd& linear)
{
  linear_ = linear;
}

void modus::GeometricJacobian::SetRotationalJacobian(const Eigen::MatrixXd& rotational)
{
  rotational_ = rotational;
}

Eigen::RowVectorXd modus::GeometricJacobian::GetGradientAtPointInDirection
  (const Eigen::Vector3d& point, const Eigen::Vector3d& direction)
{
  // Compute the linear velocity jacobian at a point r offset from the center of
  // the frame t. It should satisfy the equation
  // 
  //    ̇r = ̇t(θ) + (r - t) × ω(θ)
  // 
  // where ω is the angular velocity.
  Transform* T = link_->Get<Transform>();
  MODUS_ASSERT(T, "Error: link is missing Transform aspect");
  Eigen::Vector3d t = T->GetTranslation();
  Eigen::Vector3d dr = point - t;
  Eigen::Matrix3d dr_hat;
  dr_hat << 0, -dr.z(), dr.y(),
            dr.z(), 0, -dr.x(),
            -dr.y(), dr.x(), 0;
  Eigen::MatrixXd J = linear_ + dr_hat * rotational_;
  return direction.transpose() * J;
}