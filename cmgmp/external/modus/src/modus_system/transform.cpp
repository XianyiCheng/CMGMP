#include <modus/system/transform.hpp>
#include <modus/common/linear_algebra.hpp>


using namespace modus;

modus::Transform::Transform() {
  rotation_.setIdentity();
  translation_.setZero();
}

void modus::Transform::SetRotation(const Eigen::Matrix3d& rotation) {
  rotation_ = rotation;
}

Eigen::Matrix3d modus::Transform::GetRotation() {
  return rotation_;
}

void modus::Transform::SetTranslation(const Eigen::Vector3d& translation) {
  translation_ = translation;
}

Eigen::Vector3d modus::Transform::GetTranslation() {
  return translation_;
}

void modus::Transform::Invert() {
  InverseTransform(rotation_, translation_);
}