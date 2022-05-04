#include <modus/kinematics/jacobian_floating.hpp>
#include <modus/kinematics/jacobian_body.hpp>
#include <modus/system/transform.hpp>

using namespace modus;


modus:: FloatingDofsJacobian:: FloatingDofsJacobian(Link* link)
{
  link_ = link;
}

Eigen::RowVectorXd modus:: FloatingDofsJacobian::GetGradientAtPointInDirection
  (const Eigen::Vector3d& point, const Eigen::Vector3d& direction)
{
  // Get transform of link.
  Transform* T = link_->Get<Transform>();
  MODUS_ASSERT(T, "Error: link is missing Transform aspect");

  // Use body twist as the Dofs.
  BodyJacobian J(link_);
  J.SetBodyJacobian(Eigen::MatrixXd::Identity(6,6));
  return J.GetGradientAtPointInDirection(point, direction);
}