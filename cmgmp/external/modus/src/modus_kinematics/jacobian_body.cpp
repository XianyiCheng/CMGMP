#include <modus/kinematics/jacobian_body.hpp>
#include <modus/system/transform.hpp>
#include <modus/common/linear_algebra.hpp>

using namespace modus;


modus::BodyJacobian:: BodyJacobian(Link* link)
{
  link_ = link;
}

Eigen::RowVectorXd modus::BodyJacobian::GetGradientAtPointInDirection
  (const Eigen::Vector3d& point, const Eigen::Vector3d& direction)
{
  // Get transform of link.
  Transform* T = link_->Get<Transform>();
  MODUS_ASSERT(T, "Error: link is missing Transform aspect");

  // Convert the body jacobian into a distance gradient at the input point along
  // the input direction. First, transform the point and direction into the body
  // frame. Then use the twist to point velocity equation ω × p + v.
  Eigen::Matrix3d R_inv = T->GetRotation();
  Eigen::Vector3d t_inv = T->GetTranslation();
  InverseTransform(R_inv, t_inv);
  Eigen::Vector3d point_body = R_inv * point + t_inv;
  Eigen::Vector3d direction_body = R_inv * direction;
  Eigen::Matrix3d P = Adjoint(point_body);
  Eigen::MatrixXd V = body_.topRows(3) + P.transpose() * body_.bottomRows(3);
  return direction_body.transpose() * V;
}