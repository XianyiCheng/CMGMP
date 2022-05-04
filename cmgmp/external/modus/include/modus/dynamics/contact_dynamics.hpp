#pragma once
#include <modus/common/eigen.hpp>
#include <modus/system/system.hpp>
#include <modus/collision/collision.hpp>
#include <modus/kinematics/kinematics.hpp>
#include <modus/system/aspect.hpp>
#include <modus/dynamics/coefficient_of_friction.hpp>

namespace modus
{

class ContactDynamics : public Aspect {
 protected:
  System*     system_;
  State*      state_;
  Collision*  collision_;
  Kinematics* kinematics_;

 public:
  ContactDynamics(System* system);
  // ContactDynamics(Collision* collision, Kinematics* kinematics);

  size_t NumContactForces();

  // Create polyhedral contact force generators. The generators are the ordered
  // counter-clockwise in the contact frame, i.e. they are +0, 0+, -0, 0-.
  Eigen::MatrixXd GetContactForceGeneratorMatrix();

  // Get upper bounds for the contact forces given the input cs mode and ss
  // mode. The lower bounds are always zero.
  Eigen::VectorXd GetContactForceUpperBounds(const std::string& cs_mode, 
                                             const std::string& ss_mode);

  // Create unit contact force directions in the contact frame for +0, 0+, -0,
  // 0- using the average coefficient of friction between the two contacting
  // links.
  Eigen::MatrixXd CreateContactForceRays(ContactPtr contact);

  // Create contact force generator for given force direction.
  Eigen::VectorXd CreateContactForceGenerator(ContactPtr contact, 
                                              const Eigen::Vector3d& ray);
};
MODUS_DEFINE_SHARED(ContactDynamics);

}