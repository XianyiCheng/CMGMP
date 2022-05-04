#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <modus/common/eigen.hpp>
#include <modus/collision/contact.hpp>
#include <modus/modes/constraints.hpp>
#include <modus/common/assert.hpp>

namespace modus
{

// struct ContactGroups {
//   std::vector<Contacts>     groups;
//   std::map<ContactPtr, int> contactToGroupId;
// };

struct ContactFrame {

};

using ContactFramePtr = std::shared_ptr<ContactFrame>;

using Jacobians = std::map<std::pair<int,int>, Eigen::MatrixXd>;

class ContactPreprocessor {
 public:
  Contacts      contacts_;
  ContactGroups points_;
  ContactGroups parallels_;
  ContactGroups manifolds_;
  Jacobians     jacobians_;
  double        epsilon_;

  bool verbose_;

  void SetJacobians(Jacobians jacobians) { jacobians_ = jacobians; }
  void SetContacts(const Contacts& contacts) { contacts_ = contacts; }
  void SetVerbose(bool verbose) { verbose_ = verbose; }
  void SetEpsilon(double eps) { epsilon_ = eps; }

  Contacts GetContacts();
  Eigen::MatrixXd GetNormalConstraints();
  Eigen::MatrixXd GetTangentConstraints();

  /**
   * @brief This function preprocesses the input contact points to remove
   * various sources of numerical instability. Contact point locations and
   * normals are aligned. In particular, it takes the following steps
   *
   *  1. Group nearly identical contact points (i.e. within epsilon threshold).
   *  2. Group nearly parallel contact normals and exactly align their normals.
   *  3. Group nearly coplanar contacts into manifolds and consistently orient
   *     their normals.
   *  4. Reproject nearly coplanar contacts to be exactly coplanar.
   *  5. Reproject tangent directions into the plane.
   *
   * The input epsilon threshold should be larger for this function compared to
   * the enumeration functions. A reasonable value is larger than the
   * penetration error in the collision checking system. For bullet, we use
   * 1e-4.
   * 
   */
  void PreprocessContactPoints();

  /**
   * @brief This function preprocesses the contact frames into contact
   * constraints. 
   *
   */
  void PreprocessContactFrames();

  // Helper functions for preprocessing contact normals.
  void GroupNearlyParallelContactNormals();
  void AlignNearlyParallelContactNormals();

  // Helper functions for preprocessing contact manifolds.
  void GroupNearlyIdenticalContactPoints();
  void GroupNearlyCoplanarContactManifoldsAndReorient();
  void ReprojectNearlyCoplanarContactManifolds();
  void ReprojectContactPointGroup(Contacts point_group);
  void ReprojectTangentDirections();

  // Helper functions for preprocessing contact frames.
  void GroupCoplanarNormalsAndTangents();
};

// class MatroidPreprocessor {
//  public:
//   void ReduceConstraints(ContactConstraintsPtr constraints);
// };

}