#include <modus/kinematics/jacobian_spatial.hpp>
#include <modus/common/linear_algebra.hpp>

using namespace modus;


Eigen::RowVectorXd modus::SpatialJacobian::GetGradientAtPointInDirection
  (const Eigen::Vector3d& point, const Eigen::Vector3d& direction)
{
  // This function converts the spatial jacobian into a distance gradient at the
  // input point along the input direction. We do so using the twist to point
  // velocity equation ω × p + v.
  Eigen::Matrix3d P = Adjoint(point);
  Eigen::MatrixXd V = spatial_.topRows(3) + P.transpose() * spatial_.bottomRows(3);
  return direction.transpose() * V;
}