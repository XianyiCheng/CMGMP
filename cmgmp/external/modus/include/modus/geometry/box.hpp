#include <modus/common/eigen.hpp>
#include <modus/geometry/halfedgemesh.hpp>

namespace modus
{

class BoxGeometry : public HalfedgeMeshGeometry {
 protected:
  double x_, y_, z_;

 public:
  BoxGeometry(double x, double y, double z, Link* link);

  int GetType() override { return BOX; }
};
MODUS_DEFINE_SHARED(BoxGeometry);

}