#include <modus/modes/preprocess.hpp>
#include <modus/common/linear_algebra.hpp>


void modus::Contact::SwapAB() {
  std::swap(body_A, body_B);
  std::swap(link_A, link_B);
  point_A.swap(point_B);
  contact_normal_B *= -1;
  tangent_x.swap(tangent_y);
}

void modus::ContactPreprocessor::GroupNearlyIdenticalContactPoints()
{
  // Groups points that are within epsilon of each other.
  points_ = {{contacts_[0]}};
  for (size_t i = 1; i < contacts_.size(); i++) {
    ContactPtr c = contacts_[i];
    bool matched = false;
    for (size_t j = 0; j < points_.size(); j++) {
      ContactPtr p = points_[j][0];
      if (((c->point_A - p->point_A).norm() < epsilon_) ||
          ((c->point_B - p->point_A).norm() < epsilon_) ||
          ((c->point_A - p->point_B).norm() < epsilon_) ||
          ((c->point_B - p->point_B).norm() < epsilon_)) {
        c->point_id_ = j;
        points_[j].push_back(c);
        matched = true;
        break;
      }
    }
    if (!matched) {
      c->point_id_ = points_.size();
      points_.push_back({c});
    }
  }
  if (verbose_) {
    for (size_t i = 0; i < points_.size(); i++) {
      if (points_[i].size() > 1) {
        std::cout << "Identified contact point..." << std::endl;
        for (size_t j = 0; j < points_[i].size(); j++) {
          std::cout << *points_[i][j] << std::endl;
        }
      }
    }
  }
}

void modus::ContactPreprocessor::GroupNearlyParallelContactNormals()
{
  // Group the parallel manifolds. Note, normal vector orientation does not
  // matter.
  parallels_ = {{contacts_[0]}};
  for (size_t i = 1; i < contacts_.size(); i++) {
    ContactPtr c = contacts_[i];
    bool shares_normal = false;
    for (size_t j = 0; j < parallels_.size(); j++) {
      Eigen::Vector3d n1 = parallels_[j][0]->contact_normal_B;
      Eigen::Vector3d n2 = c->contact_normal_B;
      n1.normalize();
      n2.normalize();
      if (n1.dot(n2) < 0) {
        n2 *= -1;
      }
      if ((n1 - n2).norm() < epsilon_) {
        parallels_[j].push_back(c);
        shares_normal = true;
        break;
      }
    }
    if (!shares_normal) {
      parallels_.push_back({c});
    }
  }
  {
    size_t n = 0;
    for (size_t i = 0; i < parallels_.size(); i++) {
      n += parallels_[i].size();
    }
    MODUS_ASSERT(n == contacts_.size());
  }
}

void modus::ContactPreprocessor::AlignNearlyParallelContactNormals()
{
  // Align parallel contact manifolds.
  for (size_t i = 0; i < parallels_.size(); i++) {
    // Average normal vector.
    Eigen::Vector3d n = Eigen::Vector3d::Zero();
    // Reference normal vector.
    Eigen::Vector3d n0 = parallels_[i][0]->contact_normal_B;
    for (size_t j = 0; j < parallels_[i].size(); j++) {
      // Get j-th normal vector.
      const Eigen::Vector3d& n_j = parallels_[i][j]->contact_normal_B;
      if (n0.dot(n_j) < 0) {
        n += -n_j;
      } else {
        n += n_j;
      }
    }
    // Take average.
    n /= (double) parallels_[i].size();
    // Set normal vectors to average.
    for (size_t j = 0; j < parallels_[i].size(); j++) {
      const Eigen::Vector3d& n_j = parallels_[i][j]->contact_normal_B;
      if (n0.dot(n_j) < 0) {
        parallels_[i][j]->contact_normal_B = -n;
      } else {
        parallels_[i][j]->contact_normal_B =  n;
      }
    }
  }
}

void modus::ContactPreprocessor::GroupNearlyCoplanarContactManifoldsAndReorient()
{
  // Step through each set of contacts and group by dot product
  // distance.
  manifolds_ = {{contacts_[0]}};
  for (size_t i = 1; i < contacts_.size(); i++) {
    ContactPtr c = contacts_[i];
    bool in_a_manifold = false;
    for (size_t j = 0; j < manifolds_.size(); j++) {
      // Get manifold's A and B distance along the manifold's normal.
      ContactPtr c_m = manifolds_[j][0];
      const Eigen::Vector3d& n_m = c_m->contact_normal_B;
      MODUS_ASSERT(std::abs(n_m.norm() - 1) < 1e-8);
      double d_A_m = c_m->point_A.dot(n_m);
      double d_B_m = c_m->point_B.dot(n_m);
      MODUS_ASSERT(std::abs(d_A_m - d_B_m) < epsilon_,
        "Warning: input contact point has separation distance %f > %f", 
        c_m->contact_distance, epsilon_);

      // Check if the current normal is aligned with the manifold.
      const Eigen::Vector3d& n_c = c->contact_normal_B;
      MODUS_ASSERT(std::abs(n_c.norm() - 1) < 1e-8);
      if (n_c.dot(n_m) < 0) {
        // This should reorient n_c. But we test just to be sure.
        c->SwapAB();
        MODUS_ASSERT(n_c.dot(n_m) > 0);
      }
      // This value should actually be 0 if they are on the same manifold due to
      // the previous alignment of nearly parallel normals.
      if ((n_c - n_m).norm() > epsilon_) {
        continue;
      }

      // Check if the current contact point is coplanar with the manifold
      // points.
      double d_A = c->point_A.dot(n_m);
      double d_B = c->point_B.dot(n_m);
      // TODO FIXME should this be && instead?
      if (std::abs(d_A - d_A_m) > epsilon_ ||
          std::abs(d_B - d_B_m) > epsilon_) {
        continue;
      }

      // Finally, add the contact to the manifold.
      manifolds_[j].push_back(c);
      in_a_manifold = true;
      break;
    }
    if (!in_a_manifold) {
      manifolds_.push_back({c});
    }
  }
  if (verbose_) {
    for (size_t i = 0; i < manifolds_.size(); i++) {
      std::cout << "Manifold..." << std::endl;
      for (size_t j = 0; j < manifolds_[i].size(); j++) {
        std::cout << *manifolds_[i][j] << std::endl;
      }
    }
  }
}

void modus::ContactPreprocessor::ReprojectNearlyCoplanarContactManifolds()
{
  // This function requires both points_ and manifolds_ be constructed.
  MODUS_ASSERT(points_.size() > 0);
  MODUS_ASSERT(manifolds_.size() > 0);
  // For each manifold, calculate the target hyperplane offset.
  for (size_t i = 0; i < manifolds_.size(); i++) {
    // Average contact point offsets along hyperplane normal.
    double avg_offset = 0;
    const Eigen::Vector3d& n = manifolds_[i][0]->contact_normal_B;
    for (size_t j = 0; j < manifolds_[i].size(); j++) {
      avg_offset += manifolds_[i][j]->point_A.dot(n);
      avg_offset += manifolds_[i][j]->point_B.dot(n);
    }
    avg_offset /= manifolds_[i].size() * 2;
    // Write the average offset to each point in the manifold.
    for (size_t j = 0; j < manifolds_[i].size(); j++) {
      manifolds_[i][j]->target_offset_ = avg_offset;
    }
  }
  // Reproject each contact point group using the average offsets.
  for (size_t i = 0; i < points_.size(); i++) {
    ReprojectContactPointGroup(points_[i]);
  }
}

void modus::ContactPreprocessor::ReprojectContactPointGroup(Contacts point_group)
{
  MODUS_ASSERT(1 <= point_group.size() && point_group.size() <= 3);
  // Solve an equality constraint quadratic program to determine the new contact
  // point location. When the # points <= 2, we can use the KKT conditions to
  // solve this. When # points == 3, we solve a linear system.
  Eigen::Vector3d x;
  size_t n_p = point_group.size();
  if (n_p <= 2) {
    // Build the KKT system, [I N^T; N 0][x; λ] = [x0; d].
    Eigen::MatrixXd N(n_p, 3);
    for (size_t i = 0; i < n_p; i++) {
      N.row(i) = point_group[i]->contact_normal_B;
    }
    Eigen::MatrixXd x0 = point_group[0]->point_A;
    Eigen::MatrixXd A(3 + n_p, 3 + n_p);
    Eigen::VectorXd b(3 + n_p);
    // A.setIdentity();
    A.setZero();
    A.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
    A.block(3, 0, n_p, 3) = N;
    A.block(0, 3, 3, n_p) = N.transpose();
    b.block(0,0,3,1) = point_group[0]->point_A;
    for (size_t i = 0; i < n_p; i++) {
      b[3 + i] = point_group[i]->target_offset_;
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
    qr.compute(A);
    Eigen::VectorXd x_lambda = qr.solve(b);
    x = x_lambda.block<3,1>(0,0);
  } else {
    Eigen::MatrixXd N(3, 3);
    for (size_t i = 0; i < 3; i++) {
      N.row(i) = point_group[i]->contact_normal_B;
    }
    Eigen::Vector3d b;
    for (size_t i = 0; i < n_p; i++) {
      b[i] = point_group[i]->target_offset_;
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
    qr.compute(N);
    x = qr.solve(b);
  }
  // Assert that the new point is at the correct target offset.
  for (size_t i = 0; i < point_group.size(); i++) {
    const Eigen::Vector3d& n = point_group[i]->contact_normal_B;
    double target = point_group[i]->target_offset_;
    MODUS_ASSERT(std::abs(x.dot(n) - target) < 1e-9);
  }
  // Assign new contact point to group.
  for (size_t i = 0; i < n_p; i++) {
    point_group[i]->point_A = x;
    point_group[i]->point_B = x;
  }
}

void modus::ContactPreprocessor::ReprojectTangentDirections()
{
  // In previous adjustments, the normal directions have probably shifted. This
  // function reprojects the tangent directions to lie in the normal plane.
  for (size_t i = 0; i < contacts_.size(); i++) {
    ContactPtr c = contacts_[i];
    const Eigen::Vector3d& n = c->contact_normal_B;
    c->tangent_x = c->tangent_y.cross(n).normalized();
    c->tangent_y = n.cross(c->tangent_x).normalized();
    MODUS_ASSERT(std::abs(c->contact_normal_B.norm() - 1) < 1e-9);
    MODUS_ASSERT(std::abs(c->tangent_x.norm() - 1) < 1e-9,
      "Error: tangent x has norm %f", c->tangent_x.norm());
    MODUS_ASSERT(std::abs(c->tangent_y.norm() - 1) < 1e-9,
      "Error: tangent y has norm %f", c->tangent_y.norm());
  }
}

void modus::ContactPreprocessor::PreprocessContactPoints()
{
  MODUS_ASSERT(contacts_.size() > 0);
  MODUS_ASSERT(epsilon_ > 0);

  GroupNearlyParallelContactNormals();
  AlignNearlyParallelContactNormals();

  GroupNearlyIdenticalContactPoints();
  GroupNearlyCoplanarContactManifoldsAndReorient();
  ReprojectNearlyCoplanarContactManifolds();
  ReprojectTangentDirections();
}
