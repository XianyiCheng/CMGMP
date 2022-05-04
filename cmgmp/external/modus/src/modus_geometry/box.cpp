#include <vector>
#include <modus/geometry/box.hpp>


using namespace modus;

modus::BoxGeometry::BoxGeometry(double x0, double y0, double z0, Link* link)
  : x_(x0), y_(y0), z_(z0), HalfedgeMeshGeometry(link)
{
  // Create box vertices.
  std::vector<Eigen::Vector3d> vertices;
  std::vector<double> X = {-x0/2, x0/2};
  std::vector<double> Y = {-y0/2, y0/2};
  std::vector<double> Z = {-z0/2, z0/2};
  int k = 0;
  for (auto x : X) {
    for (auto y : Y) {
      for (auto z : Z) {
        vertices.push_back(Eigen::Vector3d(x, y, z));
      }
    }
  }
  // Build from convex hull.
  BuildConvex(vertices);
}