#pragma once
#include <modus/common/memory.hpp>
#include <modus/common/eigen.hpp>
#include <modus/system/aspect.hpp>


namespace modus
{

class Transform : public Aspect {//, public std::enable_shared_from_this<Transform> {
 public:
  Eigen::Matrix3d rotation_;
  Eigen::Vector3d translation_;

  Transform();

  void SetRotation(const Eigen::Matrix3d& rotation);
  Eigen::Matrix3d GetRotation();

  void SetTranslation(const Eigen::Vector3d& translation);
  Eigen::Vector3d GetTranslation();

  void Invert();
  Transform Inverted();
};
MODUS_DEFINE_SHARED(Transform);

}