#include <modus/modes/constraints.hpp>


modus::ContactConstraints::ContactConstraints
  (const Contacts& contacts, const JacobianMap& jacobians, bool verbose)
{
  // Build normal velocity and tangent velocity constraints.
  size_t n_contacts = contacts.size();
  size_t n_dofs = jacobians.begin()->second->GetTotalDofs();
  N.resize(n_contacts, n_dofs);
  N.setZero();
  T.resize(n_contacts * 2, n_dofs);
  T.setZero();
  d.resize(n_contacts);
  for (size_t i = 0; i < contacts.size(); i++) {
    if (verbose) {
      std::cout << "Adding contact" << std::endl;
      std::cout << *contacts[i] << std::endl;
    }

    ContactPtr c = contacts[i];
    // Get the jacobian of the chain from body A -> link A -> point A and create
    // the normal velocity constraint.
    auto pair_A = std::make_pair(c->body_A, c->link_A);
    JacobianPtr J_A = jacobians.at(pair_A);
    if (J_A) {
      N.row(i) += J_A->GetGradientAtPointInDirectionFullDofs(c->point_A, c->contact_normal_B);
      T.row(2*i + 0) += J_A->GetGradientAtPointInDirectionFullDofs(c->point_A, c->tangent_x);
      T.row(2*i + 1) += J_A->GetGradientAtPointInDirectionFullDofs(c->point_A, c->tangent_y);
    }

    // Get the jacobian of the chain from body B -> link B -> point B and create
    // the normal velocity constraint.
    auto pair_B = std::make_pair(c->body_B, c->link_B);
    JacobianPtr J_B = jacobians.at(pair_B);
    if (J_B) {
      N.row(i) += J_B->GetGradientAtPointInDirectionFullDofs(c->point_B, c->contact_normal_B);
      T.row(2*i + 0) += J_B->GetGradientAtPointInDirectionFullDofs(c->point_B, c->tangent_x);
      T.row(2*i + 1) += J_B->GetGradientAtPointInDirectionFullDofs(c->point_B, c->tangent_y);
    }

    // Write the separating distance.
    d[i] = c->contact_distance;
  }
  N *= -1;
}

void modus::ContactConstraints::Partition(const std::string& cs_mode,
                                          const std::string& ss_mode)
{
  // Partition the normal velocity constraints based on cs mode.
  size_t n_contacting = std::count(cs_mode.begin(), cs_mode.end(), '0');
  size_t n_separating = std::count(cs_mode.begin(), cs_mode.end(), '-');
  size_t n_dofs = N.cols();
  size_t k = 0, k_eq = 0;
  MODUS_ASSERT(cs_mode.size() == n_dofs);
  N_eq.resize(n_contacting, n_dofs);
  d_eq.resize(n_contacting);
  for (size_t i = 0; i < n_dofs; i++) {
    if (cs_mode[i] == '0') {
      N_eq.row(k_eq) = N.row(i);
      d_eq[k_eq++] = d[i];
    } else if (cs_mode[i] == '-') {
      N.row(k) = N.row(i);
      d[k++] = d[i];
    } else {
      MODUS_ASSERT(false, "Error: invalid mode string %s", cs_mode.c_str());
    }
  }
  N.conservativeResize(n_separating, n_dofs);
  MODUS_ASSERT(k == n_separating);
  MODUS_ASSERT(k_eq == n_contacting);

  // Partition the tangent velocity constraints based on cs mode and ss mode.
  size_t n_partially_sticking = std::count(ss_mode.begin(), ss_mode.end(), '0');
  size_t n_sliding = ss_mode.size() - n_partially_sticking;
  size_t j = 0;
  k = 0;
  k_eq = 0;
  for (size_t i = 0; i < cs_mode.size(); i++) {
    // if cs_mode[i] == '0', add the tangent hyperplanes based on their modes
    if (cs_mode[i] == '0') {
      char m = ss_mode[j++];
      if (m == '0') {
        T_eq.row(k_eq++) = T.row(2*i+0);
      } else if (m == '+') {
        T.row(k++) = T.row(2*i+0);
      } else if (m == '-') {
        T.row(k++) = T.row(2*i+0);
      } else {
        MODUS_ASSERT(false);
      }
      m = ss_mode[j++];
      if (m == '0') {
        T_eq.row(k_eq++) = T.row(2*i+1);
      } else if (m == '+') {
        T.row(k++) = T.row(2*i+1);
      } else if (m == '-') {
        T.row(k++) = T.row(2*i+1);
      } else {
        MODUS_ASSERT(false);
      }
    }
  }
  MODUS_ASSERT(k == n_sliding);
  MODUS_ASSERT(k_eq == n_partially_sticking);
  T.conservativeResize(n_sliding, n_dofs);
}

void modus::ContactConstraintGenerator::UpdateConstraints
  (Contacts contacts, JacobianMap jacobians, double eps)
{
  dynamics_contacts_ = contacts;
  jacobians_ = jacobians;
  epsilon_ = eps;

  // Make a copy of the contact points which we will numerically stabilize by
  // shifting and rotating.
  filtered_contacts_.resize(0);
  for (size_t i = 0; i < contacts.size(); i++) {
    filtered_contacts_.push_back(ContactPtr(new Contact(*contacts[i])));
  }

  // Preprocess contact points by making them coplanar, aligning normals, and
  // rotating tangents.
  ContactPreprocessor contact_pre;
  contact_pre.SetContacts(filtered_contacts_);
  // contact_pre.SetJacobians(jacobians);
  contact_pre.SetEpsilon(epsilon_);
  contact_pre.SetVerbose(true);
  contact_pre.PreprocessContactPoints();
}

modus::ContactConstraintsPtr modus::ContactConstraintGenerator::GetConstraints
  (const std::string& cs_mode, const std::string& ss_mode, int flags)
{
  // Select appropriate level of preprocessed contacts.
  Contacts contacts;
  if (flags & CC_DYNAMICS) {
    contacts = dynamics_contacts_;
  } else if (flags & CC_CONTACT_FILTER) {
    std::cout << "do this" << std::endl;
    contacts = filtered_contacts_;
  }

  // Move contacts, normals, and tangents around to get well-conditioned
  // inputs.
  // if (flags & CC_CONTACT_FILTER) {
  //   ContactPreprocessor preprocessor;
  //   preprocessor.SetContacts(contacts);
  //   preprocessor.SetEpsilon(epsilon_);
  //   preprocessor.PreprocessContactPoints();
  // }

  // Build constraint matrices.
  ContactConstraintsPtr constraints(new ContactConstraints(contacts, jacobians_, true));

  // Use matroid projective equivalence to reduce constraint dimensions.
  if (flags & CC_MATROID_FILTER) {

  }

  return constraints;
}