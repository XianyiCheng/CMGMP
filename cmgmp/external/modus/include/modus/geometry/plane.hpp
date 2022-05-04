#pragma once
#include <modus/common/eigen.hpp>
#include <modus/geometry/halfedgemesh.hpp>


namespace modus
{

class PlaneGeometry : public HalfedgeMeshGeometry {
 protected:
  Eigen::Vector3d normal_;
  double offset_;

 public:
  PlaneGeometry(const Eigen::Vector3d& normal, double offset);

  int GetType() override { return PLANE; }

  Eigen::Vector3d GetNormal() { return normal_; }
  double GetOffset() { return offset_; }
};
MODUS_DEFINE_SHARED(PlaneGeometry);

}