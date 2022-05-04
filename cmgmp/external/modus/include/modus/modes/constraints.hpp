#pragma once
#include <modus/collision/contact.hpp>
#include <modus/kinematics/jacobian.hpp>
#include <modus/modes/preprocess.hpp>


namespace modus
{

/**
 * @brief Data structure which stores the normal and tangent velocity
 * constraints. Normal velocities constraints follow the convention Nx - d <= 0,
 * which means the signs are flipped. If Nx - d < 0, then v_p⋅n > 0.
 *
 * However, tangent velocity constraints do not flip signs. If Tx > 0, then
 * v_p⋅t_x > 0.
 *
 */
struct ContactConstraints {
  // Normal velocity constraints.
  Eigen::MatrixXd N;
  Eigen::VectorXd d;
  Eigen::MatrixXd N_eq;
  Eigen::VectorXd d_eq;
  Eigen::MatrixXd T;
  Eigen::MatrixXd T_eq;

  // Projection matrix for dimensionality reduction.
  Eigen::MatrixXd R;

  ContactConstraints(const Contacts& contacts, const JacobianMap& jacobians,
                     bool verbose=false);

  void Partition(const std::string& cs_mode,
                 const std::string& ss_mode);
};

using ContactConstraintsPtr = std::shared_ptr<ContactConstraints>;

class ContactConstraintGenerator {
 protected:
  Contacts    dynamics_contacts_;
  Contacts    filtered_contacts_;
  JacobianMap jacobians_;
  double      epsilon_;

 public:
  enum {
    // Unprocessed constraints for dynamic simulation/optimization.
    CC_DYNAMICS = 1,
    // Preprocess contacts to get well-conditioned constraints.
    CC_CONTACT_FILTER = 2,
    // Use matroid projective equivalence to reduce constraint dimensions.
    CC_MATROID_FILTER = 4,
  };

  void UpdateConstraints(Contacts contacts, JacobianMap jacobians, double eps);

  /**
   * @brief This function creates normal and tangent velocity constraint
   * matrices based on the input cs mode, ss mode, and flags.
   *
   * x = GetConstraints("", "", CC_DYNAMICS) generates N, d, T matrices for
   * optimization using rotated tangent vectors and no other changes.
   *
   * x = GetConstraints(cs_mode, ss_mode, CC_DYNAMICS) partitions N, d, T
   * matrices based on the input mode, overwriting N, d, T for the inequality
   * constraints and N_eq, d_eq, T_eq for equality constraints.
   *
   * x = GetConstraints("", "", CC_CONTACT_FILTER) uses shifted contact points
   * (coplanar, idential points associated), aligned normals, and rotated
   * tangents. Intended for enumeration only, d is always 0.
   *
   * x = GetConstraints(cs_mode, ss_mode, CC_CONTACT_FILTER) partitioned version
   * of previous call.
   *
   * x = GetConstraints("", "", CC_CONTACT_FILTER & CC_MATROID_FILTER) uses
   * preprocessed contacts but also applies matroid-based projective equivalence
   * to reduce constraint dimensionality. 
   *
   * x = GetConstraints(cs_mode, ss_mode, CC_CONTACT_FILTER & CC_MATROID_FILTER)
   * partitioned version where the matroid filter is applied after partitioning.
   *
   * @param cs_mode 
   * @param ss_mode 
   * @param flags 
   */
  ContactConstraintsPtr GetConstraints(const std::string& cs_mode, 
                                       const std::string& ss_mode, 
                                       int flags);
};

using ContactConstraintGeneratorPtr = std::shared_ptr<ContactConstraintGenerator>;

}