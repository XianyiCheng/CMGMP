#pragma once
#include <vector>
#include <modus/common/memory.hpp>
#include <modus/common/eigen.hpp>
#include <modus/system/link.hpp>
#include <modus/geometry/geometry.hpp>
#include <modus/geometry/halfedge.hpp>
#include <modus/geometry/vertex.hpp>
#include <modus/geometry/edge.hpp>
#include <modus/geometry/face.hpp>


namespace modus
{

class HalfedgeMeshGeometry : public Geometry {
 protected:
  Link* link_;

 public:
  std::vector<HalfedgePtr> halfedges_;    //
  std::vector<VertexPtr>   vertices_;     //
  std::vector<EdgePtr>     edges_;        //
  std::vector<FacePtr>     faces_;        //
  std::vector<FacePtr>     boundaries_;   //
  float margin_;

  // GLuint vao;
  // GLuint vbo;
  // GLuint ebo;
  // GLsizei num_elem_draw;

  HalfedgeMeshGeometry(Link* link = nullptr);

  int GetType() override { return MESH; }

  Transform GetTransformWorld();

  // Get vertex positions in local coordinates.
  Eigen::MatrixXd GetVertexPositions();

  // Get vertex positions in world coordinates.
  Eigen::MatrixXd GetVertexPositionsWorld();

  HalfedgePtr NewHalfedge();
  VertexPtr NewVertex();
  EdgePtr NewEdge();
  FacePtr NewFace();
  FacePtr NewBoundary();

  void Load(const std::string& path);
  void Build(const std::vector<Eigen::VectorXi>& simplices, 
             const std::vector<Eigen::Vector3d>& positions);
  void BuildConvex(const std::vector<Eigen::Vector3d>& points);
  void Reindex();
  // void Init(const Eigen::Vector3f& color);
  // virtual void Draw(int program);
  // std::vector<mesh_geom::ShapePtr> Primitives();

  // Mesh calculations.
  double SurfaceArea();
  double BarycentricDualArea(VertexPtr v);
  double AngleDefect(VertexPtr v);
  double ScalarGaussCurvature(VertexPtr v);
  double ScalarMeanCurvature(VertexPtr v);

  // Mesh processing.
  void SplitEdge(EdgePtr e);
  void FlipEdge(EdgePtr e);
  void SplitPolygon(FacePtr f);
  void Triangulate();
  void LoopSubdivide();
  void Subdivide41();
};
MODUS_DEFINE_SHARED(HalfedgeMeshGeometry);

}