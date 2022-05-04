#include <iostream>
#include <map>
#include <set>
#include <modus/geometry/halfedgemesh.hpp>
#include <modus/common/linear_algebra.hpp>
#include <modus/common/convex_hull.hpp>
// #include <modus/geometry/opengl.hpp>


using namespace modus;

HalfedgeMeshGeometry::HalfedgeMeshGeometry(Link* link) {
  link_ = link;
}

Transform HalfedgeMeshGeometry::GetTransformWorld() {
  return *link_->Get<Transform>();
}

// Get vertex positions in local coordinates.
Eigen::MatrixXd HalfedgeMeshGeometry::GetVertexPositions() {
  // Get 
  Eigen::MatrixXd positions(3, vertices_.size());
  for (size_t i = 0; i < vertices_.size(); i++) {
  positions.col(i) = vertices_[i]->position;
  }
  return positions;
}

// Get vertex positions in world coordinates.
Eigen::MatrixXd HalfedgeMeshGeometry::GetVertexPositionsWorld() {
  Transform tf = GetTransformWorld();
  Eigen::Matrix3d R = tf.GetRotation();
  Eigen::Vector3d t = tf.GetTranslation();
  Eigen::MatrixXd positions(3, vertices_.size());
  for (size_t i = 0; i < vertices_.size(); i++) {
  positions.col(i) = R * vertices_[i]->position + t;
  }
  return positions;
}

HalfedgePtr HalfedgeMeshGeometry::NewHalfedge() {
  HalfedgePtr h = std::make_shared<Halfedge>();
  halfedges_.push_back(h);
  return h;
}

VertexPtr HalfedgeMeshGeometry::NewVertex() {
  return *vertices_.insert(vertices_.end(), std::make_shared<Vertex>());
}

EdgePtr HalfedgeMeshGeometry::NewEdge() {
  return *edges_.insert(edges_.end(), std::make_shared<Edge>());
}

FacePtr HalfedgeMeshGeometry::NewFace() {
  return *faces_.insert(faces_.end(), std::make_shared<Face>());
}

FacePtr HalfedgeMeshGeometry::NewBoundary() {
  return *boundaries_.insert(boundaries_.end(), std::make_shared<Face>(true));
}

void HalfedgeMeshGeometry::Build(const std::vector<Eigen::VectorXi>& simplices, 
                 const std::vector<Eigen::Vector3d>& positions) {
  // Clear elements.
  halfedges_.clear();
  vertices_.clear();
  edges_.clear();
  faces_.clear();
  boundaries_.clear();

  // Maps a vertex index to the corresponding vertex.
  std::map<int, VertexPtr> indexToVertex;

  // Also store the vertex degree, i.e., the number of polygons that use
  // each vertex this information will be used to check that the mesh is
  // manifold.
  std::map<VertexPtr, int> vertexDegree;

  // Do some basic sanity checks on the input.
  for (Eigen::VectorXi p : simplices) {
    if (p.size() < 3) {
      throw std::runtime_error("Error converting polygons to halfedge mesh: each polygon must have at least three vertices.");
    }

    // We want to count the number of distinct vertex indices in this
    // polygon, to make sure it's the same as the number of vertices
    // in the polygon---if they disagree, then the polygon is not valid
    // (or at least, for simplicity we don't handle polygons of this type!).
    std::set<int> polygonIndices;

    // Loop over polygon vertices.
    for (int i = 0; i < p.size(); i++) {
      polygonIndices.insert(p[i]);

      // Allocate one vertex for each new index we encounter.
      if (indexToVertex.find(p[i]) == indexToVertex.end()) {
        VertexPtr v = NewVertex();
        v->halfedge = nullptr;
        indexToVertex[p[i]] = v;
        vertexDegree[v] = 1; // We've now seen this vertex only once.
      } else {
        // keep track of the number of times we've seen this vertex
        vertexDegree[indexToVertex[p[i]]]++;
      }
    } // End loop over polygon vertices.

    // Check that all the vertices of the current polygon are distinct.
    int degree = p.size();
    if (polygonIndices.size() < degree) {
      throw std::runtime_error("Error converting polygons to halfedge mesh: one of the input polygons does not have distinct vertices");
    }
  } // end basic sanity check on input

  // The number of vertices in the mesh is the number of unique indices
  // seen in the input.
  int nVertices = indexToVertex.size();

  // The number of faces is the number of polygons in the input.
  int nFaces = simplices.size();
  // Allocate storage for faces in our new mesh
  faces_.resize(nFaces);  
  for (int i = 0; i < nFaces; i++) {
    faces_[i] = std::make_shared<Face>();
  }

  // We will store a map from ordered pairs of vertex indices to
  // the corresponding halfedge object in our new (halfedge) mesh
  // this map gets constructed during the next loop over polygons.
  std::map<std::pair<int, int>, HalfedgePtr> pairToHalfedge;

  // Next, we build the halfedge connectivity by again looping over
  // polygons.
  for (int i = 0; i < nFaces; i++) {
    const Eigen::VectorXi& p = simplices[i];
    FacePtr f = faces_[i];

    std::vector<HalfedgePtr> faceHalfedges; // Cyclically ordered list of 
                        // halfedges of this face.
    int degree = p.size();                  // Number of vertices in this polgyon.

    // Loop over the halfedges of this face (ordered pairs of consecutive
    // vertices).
    for (int k = 0; k < degree; k++) {
      int a = p[k];                 // current index
      int b = p[(k + 1) % degree];  // next index, in cyclic order
      std::pair<int, int> ab(a, b);
      HalfedgePtr hab;

      // check if this halfedge already exists; if so, we have a problem!
      if (pairToHalfedge.find(ab) != pairToHalfedge.end()) {
        std::cerr << "Error converting polygons to halfedge mesh: found multiple "
            "oriented edges with indices ("
          << a << ", " << b << ")." << std::endl;
        std::cerr << "This means that either (i) more than two faces contain this "
            "edge (hence the surface is nonmanifold), or"
          << std::endl;
        std::cerr << "(ii) there are exactly two faces containing this edge, but "
            "they have the same orientation (hence the surface is"
          << std::endl;
        std::cerr << "not consistently oriented." << std::endl;
        throw std::runtime_error("Error converting polygons to halfedge mesh");
      } else  // otherwise, the halfedge hasn't been allocated yet
      {
        // so, we point this vertex pair to a new halfedge
        hab = NewHalfedge();
        pairToHalfedge[ab] = hab;

        // link the new halfedge to its face
        hab->face = f;
        hab->face->halfedge = hab;

        // also link it to its starting vertex
        hab->vertex = indexToVertex[a];
        hab->vertex->halfedge = hab;

        // keep a list of halfedges in this face, so that we can later
        // link them together in a loop (via their "next" pointers)
        faceHalfedges.push_back(hab);
      }

      // Also, check if the twin of this halfedge has already been constructed
      // (during construction of a different face).  If so, link the twins
      // together and allocate their shared halfedge.  By the end of this pass
      // over polygons, the only halfedges that will not have a twin will hence
      // be those that sit along the domain boundary.
      std::pair<int, int> ba(b, a);
      std::map<std::pair<int, int>, HalfedgePtr>::iterator iba = pairToHalfedge.find(ba);
      if (iba != pairToHalfedge.end()) {
        HalfedgePtr hba = iba->second;

        // link the twins
        hab->twin = hba;
        hba->twin = hab;

        // allocate and link their edge
        EdgePtr e = NewEdge();
        hab->edge = e;
        hba->edge = e;
        e->halfedge = hab;
      } else {  // If we didn't find a twin...
        // ...mark this halfedge as being twinless by pointing
        // it to the end of the list of halfedges. If it remains
        // twinless by the end of the current loop over polygons,
        // it will be linked to a boundary face in the next pass.
        hab->twin = nullptr;
      }

    }  // end loop over the current polygon's halfedges

    // Now that all the halfedges of this face have been allocated,
    // we can link them together via their "next" pointers.
    for (int k = 0; k < degree; k++) {
      int j = (k + 1) % degree;  // index of the next halfedge, in cyclic order
      faceHalfedges[k]->next = faceHalfedges[j];
    }

  } // done building basic halfedge connectivity

  // For each vertex on the boundary, advance its halfedge pointer to one that
  // is also on the boundary.
  for (auto v : vertices_) {
    // loop over halfedges around vertex
    HalfedgePtr h = v->halfedge;
    do {
      if (h->twin == nullptr) {
        v->halfedge = h;
        break;
      }
      h = h->twin->next;
    } while (h != v->halfedge);  // end loop over halfedges around vertex

  }  // done advancing halfedge pointers for boundary vertices

  // Construct faces for each boundary component.
  for (int i = 0; i < halfedges_.size(); i++) { // loop over all halfedges
    HalfedgePtr h = halfedges_[i];
    // Any halfedge that does not yet have a twin is on the boundary of
    // the domain. If we follow the boundary around long enough we will
    // of course eventually make a closed loop we can represent this
    // boundary loop by a new face. To make clear the distinction between
    // faces and boundary loops, the boundary face will (i) have a flag
    // indicating that it is a boundary loop, and (ii) be stored in a
    // list of boundaries, rather than the usual list of faces.  The
    // reason we need the both the flag *and* the separate list is that
    // faces are often accessed in two fundamentally different ways:
    // either by (i) local traversal of the neighborhood of some mesh
    // element using the halfedge structure, or (ii) global traversal of
    // all faces (or boundary loops).
    if (h->twin == nullptr) {
      FacePtr b = NewBoundary();
      std::vector<HalfedgePtr> boundaryHalfedges; // keep a list of halfedges along
                            // the boundary, so we can link
                            // them together

      // We now need to walk around the boundary, creating new halfedges
      // and edges along the boundary loop as we go.
      HalfedgePtr i = h;
      do {
        // create a twin, which becomes a halfedge of the boundary loop
        HalfedgePtr t = NewHalfedge();
        boundaryHalfedges.push_back(
          t);  // keep a list of all boundary halfedges, in cyclic order
        i->twin = t;
        t->twin = i;
        t->face = b;
        t->vertex = i->next->vertex;

        // create the shared edge
        EdgePtr e = NewEdge();
        e->halfedge = i;
        i->edge = e;
        t->edge = e;

        // Advance i to the next halfedge along the current boundary loop
        // by walking around its target vertex and stopping as soon as we
        // find a halfedge that does not yet have a twin defined.
        i = i->next;
        while (i != h &&  // we're done if we end up back at the beginning of
                  // the loop
            i->twin != nullptr) // otherwise, we're looking for
                      // the next twinless halfedge
                      // along the loop
        {
          i = i->twin->next;
        }
      } while (i != h);

      b->halfedge = boundaryHalfedges.front();

      // The only pointers that still need to be set are the "next" pointers of
      // the twins; these we can set from the list of boundary halfedges, but we
      // must use the opposite order from the order in the list, since the
      // orientation of the boundary loop is opposite the orientation of the
      // halfedges "inside" the domain boundary.
      int degree = boundaryHalfedges.size();
      for (int p = 0; p < degree; p++) {
        int q = (p - 1 + degree) % degree;
        boundaryHalfedges[p]->next = boundaryHalfedges[q];
      }

    } // end construction of one of the boundary loops

    // Note that even though we are looping over all halfedges, we will still
    // construct the appropriate number of boundary loops (and not, say, one
    // loop per boundary halfedge).  The reason is that as we continue to
    // iterate through halfedges, we check whether their twin has been assigned,
    // and since new twins may have been assigned earlier in this loop, we will
    // end up skipping many subsequent halfedges.

  } // done adding "virtual" faces corresponding to boundary loops

  // To make later traversal of the mesh easier, we will now advance the
  // halfedge associated with each vertex such that it refers to the *first*
  // non-boundary halfedge, rather than the last one.
  for (auto v : vertices_) {
    v->halfedge = v->halfedge->twin->next;
  }

  // Finally, we check that all vertices are manifold.
  for (auto v : vertices_) {
    // First check that this vertex is not a "floating" vertex;
    // if it is then we do not have a valid 2-manifold surface.
    if (v->halfedge == nullptr) {
      throw std::runtime_error("Error converting polygons to halfedge mesh: some vertices are "
                   "not referenced by any polygon.");
    }

    // Next, check that the number of halfedges emanating from this vertex in
    // our half edge data structure equals the number of polygons containing
    // this vertex, which we counted during our first pass over the mesh.  If
    // not, then our vertex is not a "fan" of polygons, but instead has some
    // other (nonmanifold) structure.
    int count = 0;
    HalfedgePtr h = v->halfedge;
    do {
      if (!h->face->isBoundary()) {
        count++;
      }
      h = h->twin->next;
    } while (h != v->halfedge);

    if (count != vertexDegree[v]) {
      throw std::runtime_error("Error converting polygons to halfedge mesh: at least one of the "
                   "vertices is nonmanifold.");
    }
  }  // end loop over vertices

  // Now that we have the connectivity, we copy the list of vertex
  // positions into member variables of the individual vertices.
  MODUS_ASSERT(vertices_.size() == positions.size());
  for (auto e : indexToVertex) {
    VertexPtr v = e.second;
    v->position = positions[e.first];
  }

  // Triangulate the mesh.
  Triangulate();

  // Finally, index and cache mesh elements.
  Reindex();
}

void HalfedgeMeshGeometry::BuildConvex(const std::vector<Eigen::Vector3d>& points) {
  // Build mesh from convex hull of points.
  Build(ConvexHull(points, 1e-9), points);
}

void HalfedgeMeshGeometry::Reindex() {
  for (int i = 0; i < halfedges_.size(); i++) {
    halfedges_[i]->index = i;
  }
  for (int i = 0; i < vertices_.size(); i++) {
    vertices_[i]->index = i;
  }
  for (int i = 0; i < edges_.size(); i++) {
    edges_[i]->index = i;
  }
  for (int i = 0; i < faces_.size(); i++) {
    faces_[i]->index = i;
  }
  for (int i = 0; i < boundaries_.size(); i++) {
    boundaries_[i]->index = i;
  }
}

// void HalfedgeMeshGeometry::init(const Eigen::Vector3f& color) {
//     // Get per-face vertex positions, normals, and colors.
//     int num_faces = faces_.size();
//     int num_verts = 3*num_faces;
//     int num_floats = 3*3;
//     Eigen::MatrixXf data(num_floats, num_verts);
//     int k = 0;
//     for (int i = 0; i < num_faces; i++) {
//         FacePtr f = faces_[i];
//         HalfedgePtr h = f->halfedge;
//         const Eigen::Vector3f normal = f->normal();
//         while (true) {
//             data(0,k) = h->vertex->position[0];
//             data(1,k) = h->vertex->position[1];
//             data(2,k) = h->vertex->position[2];
//             data(3,k) = normal[0];
//             data(4,k) = normal[1];
//             data(5,k) = normal[2];
//             data(6,k) = color[0];
//             data(7,k) = color[1];
//             data(8,k) = color[2];
//             h = h->next;
//             k += 1;
//             if (h == f->halfedge) {
//                 break;
//             }
//         }
//     }
//     // num_elem_draw = data.size();
//     num_elem_draw = data.cols();

//     // Setup OpenGL VAO.
//     glGenVertexArrays(1, &this->vao);
//     glGenBuffers(1, &vbo);
  
//     glBindVertexArray(vao);

//     glBindBuffer(GL_ARRAY_BUFFER, vbo);
//     glBufferData(GL_ARRAY_BUFFER, data.size()*sizeof(float), data.data(), GL_STATIC_DRAW);

//     // vertex
//     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, num_floats*sizeof(float), (void*)0);
//     glEnableVertexAttribArray(0);
//     // normal
//     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, num_floats*sizeof(float), (void*)12);
//     glEnableVertexAttribArray(1);
//     // color
//     glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, num_floats*sizeof(float), (void*)24);
//     glEnableVertexAttribArray(2);

//     glBindVertexArray(0);
// }

// void HalfedgeMeshGeometry::draw(int program) {
//     glBindVertexArray(vao);
//     glDrawArrays(GL_TRIANGLES, 0, num_elem_draw);
//     glBindVertexArray(0);
// }

// std::vector<modus::ShapePtr> HalfedgeMeshGeometry::primitives() {
//     std::vector<modus::ShapePtr> prims(faces.begin(), faces.end());
//     return prims;
// }

void HalfedgeMeshGeometry::SplitPolygon(FacePtr f) {
  // Triangulate a polygonal face.

  // Collect mesh elements.
  std::vector<HalfedgePtr> h_int;
  HalfedgePtr h0 = f->halfedge;
  HalfedgePtr h = h0;
  while (true) {
    h_int.push_back(h);
    h = h->next;
    if (h == h0)
      break;
  }
  int n = h_int.size();
  if (n == 3)
    return;
  HalfedgePtr h1 = h0->next;
  HalfedgePtr h2 = h1->next;
  HalfedgePtr h3 = h_int[n-1];
  EdgePtr e0 = h0->edge;
  EdgePtr e1 = h1->edge;
  EdgePtr e2 = h2->edge;
  EdgePtr e3 = h3->edge;
  VertexPtr v0 = h0->vertex;
  VertexPtr v1 = h1->vertex;
  VertexPtr v2 = h2->vertex;

  // Allocate mesh elements.
  FacePtr f1 = NewFace();
  HalfedgePtr h4 = NewHalfedge();
  HalfedgePtr h5 = NewHalfedge();
  EdgePtr e4 = NewEdge();

  // Assign mesh elements.
  h0->face = f1;
  h1->face = f1;
  h1->next = h4;
  h4->setNeighbors(h0, h5, v2, e4, f1);
  e4->halfedge = h4;
  h5->setNeighbors(h2, h4, v0, e4, f);
  f1->halfedge = h0;
  h3->next = h5;
  f->halfedge = h3;

  // Recurse.
  SplitPolygon(f);
}

void HalfedgeMeshGeometry::Triangulate() {
  int n_faces = faces_.size();
  for (int i = 0; i < n_faces; i++) {
    FacePtr f = faces_[i];
    SplitPolygon(f);
  }
}

