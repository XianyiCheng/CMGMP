#include <modus/geometry/plane.hpp>

using namespace modus;

modus::PlaneGeometry::PlaneGeometry(const Eigen::Vector3d& normal,double offset)
  : normal_(normal), offset_(offset) 
{
}

