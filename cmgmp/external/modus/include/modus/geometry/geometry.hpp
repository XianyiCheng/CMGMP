#pragma once
#include <modus/common/memory.hpp>
#include <modus/common/eigen.hpp>
#include <modus/system/aspect.hpp>
#include <modus/system/transform.hpp>
#include <modus/geometry/bounds.hpp>


namespace modus {

enum GeometryType {
  BOX       = 1,
  CYLINDER  = 2,
  CAPSULE   = 3,
  SPHERE    = 4,
  PLANE     = 5,
  MESH      = 6
};

class Geometry : public Aspect {
 public:

  virtual int GetType() = 0;

  virtual Transform GetTransformWorld() = 0;

  // virtual Bounds3 worldBound() = 0;
  // virtual float distance2(const Eigen::Vector3f& point) = 0;
  // virtual Eigen::Vector3f position() = 0;
};
MODUS_DEFINE_SHARED(Geometry);

}