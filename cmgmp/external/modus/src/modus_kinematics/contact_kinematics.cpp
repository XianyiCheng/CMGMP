#include <modus/kinematics/contact_kinematics.hpp>
#include <modus/kinematics/jacobian.hpp>
#include <modus/collision/collision.hpp>


using namespace modus;

modus::ContactKinematics::ContactKinematics(System* system)
  : system_(system)
{
}

size_t modus::ContactKinematics::NumNormalVelocityConstraints()
{
  Collision* collision = system_->Get<Collision>();
  MODUS_ASSERT(collision, "Error: system is missing Collision aspect");
  return collision->GetContactPoints().size();
}

size_t modus::ContactKinematics::NumTangentVelocityConstraints()
{
  Collision* collision = system_->Get<Collision>();
  MODUS_ASSERT(collision, "Error: system is missing Collision aspect");
  return 2 * collision->GetContactPoints().size();
}

ContactConstraints modus::ContactKinematics::GetNormalVelocityConstraints
  (const std::string& cs_mode, bool filter_contacts)
{
  // Get contact points.
  Collision* collision = system_->Get<Collision>();
  MODUS_ASSERT(collision, "Error: system is missing Collision aspect");
  Contacts contacts = collision->GetContactPoints();

  // Get number of dofs.
  size_t n_dofs = system_->Get<State>()->NumDofs();

  // Create normal velocity constraints.
  ContactConstraints constraints;
  for (ContactPtr c : contacts) {
    ContactConstraintPtr h(new ContactConstraint);
    // Get jacobians.
    Jacobian* J_A = system_->GetBody(c->body_A)->GetLink(c->link_A)->
                      Get<Jacobian>();
    Jacobian* J_B = system_->GetBody(c->body_B)->GetLink(c->link_B)->
                      Get<Jacobian>();
    MODUS_ASSERT(J_A || J_B);
    // Resize constraint hyperplane normal vector.
    h->n.setZero(n_dofs);
    // Create normal velocity constraint for body A.
    if (J_A) {
      Body* body = system_->GetBody(c->body_A);
      State* state = body->Get<State>();
      size_t n_body_dofs = state->NumDofs();
      size_t first_index = state->GetFirstDofIndex();
      Eigen::Ref<Eigen::RowVectorXd> qd = h->n.block(0, first_index,
                                                     1, n_body_dofs);
      qd += J_A->GetGradientAtPointInDirection(c->point_A, 
                                               c->contact_normal_B);
    }
    // Create normal velocity constraint for body B.
    if (J_B) {
      Body* body = system_->GetBody(c->body_B);
      State* state = body->Get<State>();
      size_t n_body_dofs = state->NumDofs();
      size_t first_index = state->GetFirstDofIndex();
      Eigen::Ref<Eigen::RowVectorXd> qd = h->n.block(0, first_index,
                                                     1, n_body_dofs);
      qd += J_B->GetGradientAtPointInDirection(c->point_B, 
                                               c->contact_normal_B);
    }
    h->lb = -c->contact_distance;
    h->ub = std::numeric_limits<double>::infinity();
    h->source = c;
    h->type = ContactConstraint::CC_NORMAL;
    // Add constraint.
    constraints.push_back(h);
  }
  MODUS_ASSERT(cs_mode.size() == 0 || constraints.size() == cs_mode.size());

  // Restrict constraints to the input cs mode.
  for (int i = 0; i < cs_mode.size(); i++) {
    if (cs_mode[i] == '0') {
      constraints[i]->ub = constraints[i]->lb;
    } else if (cs_mode[i] == '-') {
      // Do nothing.
    } else if (cs_mode[i] == '+') {
      MODUS_ASSERT(false);
    } else {
      MODUS_ASSERT(false);
    }
  }

  return constraints;
}

ContactConstraints modus::ContactKinematics::GetTangentVelocityConstraints
  (const std::string& cs_mode, const std::string& ss_mode, 
   bool filtered_contacts)
{
  // Get contact points.
  Collision* collision = system_->Get<Collision>();
  MODUS_ASSERT(collision, "Error: system is missing Collision aspect");
  Contacts contacts = collision->GetContactPoints();

  // Get number of dofs.
  size_t n_dofs = system_->Get<State>()->NumDofs();

  // Create tangent velocity constraints.
  ContactConstraints constraints;
  for (ContactPtr c : contacts) {
    ContactConstraintPtr h_x(new ContactConstraint);
    ContactConstraintPtr h_y(new ContactConstraint);
    // Get jacobians.
    Jacobian* J_A = system_->GetBody(c->body_A)->GetLink(c->link_A)->
                      Get<Jacobian>();
    Jacobian* J_B = system_->GetBody(c->body_B)->GetLink(c->link_B)->
                      Get<Jacobian>();
    MODUS_ASSERT(J_A || J_B);
    // Resize constraint hyperplane normal vector.
    h_x->n.setZero(n_dofs);
    h_y->n.setZero(n_dofs);
    // Create tangent velocity constraints for body A.
    if (J_A) {
      Body* body = system_->GetBody(c->body_A);
      State* state = body->Get<State>();
      size_t n_body_dofs = state->NumDofs();
      size_t first_index = state->GetFirstDofIndex();
      Eigen::Ref<Eigen::RowVectorXd> dofs_x_A = h_x->n.block(0, first_index,
                                                             1, n_body_dofs);
      Eigen::Ref<Eigen::RowVectorXd> dofs_y_A = h_y->n.block(0, first_index,
                                                             1, n_body_dofs);
      dofs_x_A += J_A->GetGradientAtPointInDirection(c->point_A, c->tangent_x);
      dofs_y_A += J_A->GetGradientAtPointInDirection(c->point_A, c->tangent_y);
    }
    // Create tangent velocity constraints for body B.
    if (J_B) {
      Body* body = system_->GetBody(c->body_A);
      State* state = body->Get<State>();
      size_t n_body_dofs = state->NumDofs();
      size_t first_index = state->GetFirstDofIndex();
      Eigen::Ref<Eigen::RowVectorXd> dofs_x_B = h_x->n.block(0, first_index,
                                                             1, n_body_dofs);
      Eigen::Ref<Eigen::RowVectorXd> dofs_y_B = h_y->n.block(0, first_index,
                                                             1, n_body_dofs);
      dofs_x_B += J_B->GetGradientAtPointInDirection(c->point_B, c->tangent_x);
      dofs_y_B += J_B->GetGradientAtPointInDirection(c->point_B, c->tangent_y);
    }
    // Add constraints.
    h_x->lb =-std::numeric_limits<double>::infinity();
    h_x->ub = std::numeric_limits<double>::infinity();
    h_x->source = c;
    h_x->type = ContactConstraint::CC_TANGENT_X;
    constraints.push_back(h_x);
    h_y->lb =-std::numeric_limits<double>::infinity();
    h_y->ub = std::numeric_limits<double>::infinity();
    h_y->source = c;
    h_y->type = ContactConstraint::CC_TANGENT_Y;
    constraints.push_back(h_y);
  }
  // MODUS_ASSERT(cs_mode.size() == 0 || constraints.size() == cs_mode.size());

  // Restrict constraints to the input ss mode.
  int j = 0;
  for (int i = 0; i < cs_mode.size(); i++) {
    if (cs_mode[i] == '0') {
      char m = ss_mode[j++];
      if (m == '0') {
        constraints[2*i+0]->lb = 0;
        constraints[2*i+0]->ub = 0;
      } else if (m == '-') {
        constraints[2*i+0]->ub = 0;
      } else if (m == '+') {
        constraints[2*i+0]->lb = 0;
      } else {
        MODUS_ASSERT(false);
      }
      m = ss_mode[j++];
      if (m == '0') {
        constraints[2*i+1]->lb = 0;
        constraints[2*i+1]->ub = 0;
      } else if (m == '-') {
        constraints[2*i+1]->ub = 0;
      } else if (m == '+') {
        constraints[2*i+1]->lb = 0;
      } else {
        MODUS_ASSERT(false);
      }
    }
  }

  return constraints;
}