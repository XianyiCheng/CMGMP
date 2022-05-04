#include <modus/system/system.hpp>


using namespace modus;

modus::System::System()
{
  // Give this system a state aspect.
  Create<State>(new SystemState(this));
  // Give this system an input aspect.
  Create<Input>(new SystemInput(this));
}

size_t modus::System::NumBodies() {
  return bodies_.size();
}

Body* modus::System::GetBody(int body_id) {
  MODUS_ASSERT(0 <= body_id && body_id < bodies_.size());
  return bodies_[body_id].get();
}

Body* modus::System::CreateBody() {
  std::unique_ptr<Body> body = std::make_unique<Body>();
  body->SetIndex(bodies_.size());
  Body* ret = body.get();
  bodies_.push_back(std::move(body));
  return ret;
}

modus::SystemState::SystemState(System* system) 
  : system_(system) {
}

size_t modus::SystemState::NumDofs() {
  // Compute the total number of dofs in this system.
  size_t n_dofs = 0;
  for (size_t i = 0; i < system_->NumBodies(); i++) {
    Body* body = system_->GetBody(i);
    State* state = body->Get<State>();
    if (!state) { continue; }
    n_dofs += state->NumDofs();
  }
  return n_dofs;
}

size_t modus::SystemState::GetFirstDofIndex() {
  return 0;
}

Eigen::VectorXd modus::SystemState::GetPositions() {
  // Aggregate the state of bodies in this system.
  size_t n_dofs = NumDofs();
  Eigen::VectorXd q(n_dofs);
  Eigen::VectorXi filled = Eigen::VectorXi::Zero(n_dofs);
  for (size_t i = 0; i < system_->NumBodies(); i++) {
    Body* body = system_->GetBody(i);
    State* state = body->Get<State>();
    if (!state) { continue; }
    size_t j = state->GetFirstDofIndex();
    size_t d = state->NumDofs();
    q.block(j, 0, d, 1) = state->GetPositions();
    MODUS_ASSERT(filled.block(j, 0, d, 1).sum() == 0);
    filled.block(j, 0, d, 1).setOnes();
  }
  return q;
}

void modus::SystemState::SetPositions(const Eigen::VectorXd& q) {
  // Set positions into the bodies of this sytem.
  size_t n_dofs = NumDofs();
  MODUS_ASSERT(q.size() == n_dofs);
  for (size_t i = 0; i < system_->NumBodies(); i++) {
    Body* body = system_->GetBody(i);
    State* state = body->Get<State>();
    if (!state) { continue; }
    size_t j = state->GetFirstDofIndex();
    size_t d = state->NumDofs();
    state->SetPositions(q.block(j, 0, d, 1));
  }
}

Eigen::VectorXd modus::SystemState::GetVelocities() {
  // Aggregate the state of bodies in this system.
  size_t n_dofs = NumDofs();
  Eigen::VectorXd qd(n_dofs);
  Eigen::VectorXi filled = Eigen::VectorXi::Zero(n_dofs);
  for (size_t i = 0; i < system_->NumBodies(); i++) {
    Body* body = system_->GetBody(i);
    State* state = body->Get<State>();
    if (!state) { continue; }
    size_t j = state->GetFirstDofIndex();
    size_t d = state->NumDofs();
    qd.block(j, 0, d, 1) = state->GetVelocities();
    MODUS_ASSERT(filled.block(j, 0, d, 1).sum() == 0);
    filled.block(j, 0, d, 1).setOnes();
  }
  return qd;
}

void modus::SystemState::SetVelocities(const Eigen::VectorXd& qd) {
  // Set positions into the bodies of this sytem.
  size_t n_dofs = NumDofs();
  MODUS_ASSERT(qd.size() == n_dofs);
  for (size_t i = 0; i < system_->NumBodies(); i++) {
    Body* body = system_->GetBody(i);
    State* state = body->Get<State>();
    if (!state) { continue; }
    size_t j = state->GetFirstDofIndex();
    size_t d = state->NumDofs();
    state->SetVelocities(qd.block(j, 0, d, 1));
  }
}

Eigen::VectorXd modus::SystemState::GetAccelerations() {
  // Aggregate the state of bodies in this system.
  size_t n_dofs = NumDofs();
  Eigen::VectorXd qdd(n_dofs);
  Eigen::VectorXi filled = Eigen::VectorXi::Zero(n_dofs);
  for (size_t i = 0; i < system_->NumBodies(); i++) {
    Body* body = system_->GetBody(i);
    State* state = body->Get<State>();
    if (!state) { continue; }
    size_t j = state->GetFirstDofIndex();
    size_t d = state->NumDofs();
    qdd.block(j, 0, d, 1) = state->GetAccelerations();
    MODUS_ASSERT(filled.block(j, 0, d, 1).sum() == 0);
    filled.block(j, 0, d, 1).setOnes();
  }
  return qdd;
}

void modus::SystemState::SetAccelerations(const Eigen::VectorXd& qdd) {
  // Set positions into the bodies of this sytem.
  size_t n_dofs = NumDofs();
  MODUS_ASSERT(qdd.size() == n_dofs);
  for (size_t i = 0; i < system_->NumBodies(); i++) {
    Body* body = system_->GetBody(i);
    State* state = body->Get<State>();
    if (!state) { continue; }
    size_t j = state->GetFirstDofIndex();
    size_t d = state->NumDofs();
    state->SetAccelerations(qdd.block(j, 0, d, 1));
  }
}

modus::SystemInput::SystemInput(System* system) 
  : system_(system) 
{
}

size_t modus::SystemInput::NumInputs() {
  // Compute the total number of inputs for this system.
  size_t n_inputs = 0;
  for (size_t i = 0; i < system_->NumBodies(); i++) {
    Body* body = system_->GetBody(i);
    Input* input = body->Get<Input>();
    if (!input) { continue; }
    n_inputs += input->NumInputs();
  }
  return n_inputs;
}

size_t modus::SystemInput::GetFirstInputIndex() {
  return 0;
}

Eigen::VectorXd modus::SystemInput::GetInputs() {
  // Aggregate the inputs of bodies in this system.
  size_t n_inputs = NumInputs();
  Eigen::VectorXd u(n_inputs);
  Eigen::VectorXi filled = Eigen::VectorXi::Zero(n_inputs);
  for (size_t i = 0; i < system_->NumBodies(); i++) {
    Body* body = system_->GetBody(i);
    Input* input = body->Get<Input>();
    if (!input) { continue; }
    size_t j = input->GetFirstInputIndex();
    size_t d = input->NumInputs();
    u.block(j, 0, d, 1) = input->GetInputs();
    MODUS_ASSERT(filled.block(j, 0, d, 1).sum() == 0);
    filled.block(j, 0, d, 1).setOnes();
  }
  return u;
}

void modus::SystemInput::SetInputs(const Eigen::VectorXd& u) {
  // Set inputs into the bodies of this sytem.
  size_t n_inputs = NumInputs();
  MODUS_ASSERT(u.size() == n_inputs);
  for (size_t i = 0; i < system_->NumBodies(); i++) {
    Body* body = system_->GetBody(i);
    Input* input = body->Get<Input>();
    if (!input) { continue; }
    size_t j = input->GetFirstInputIndex();
    size_t d = input->NumInputs();
    input->SetInputs(u.block(j, 0, d, 1));
  }
}