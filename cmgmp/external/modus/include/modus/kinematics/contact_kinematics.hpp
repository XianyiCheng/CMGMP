#pragma once
#include <modus/system/system.hpp>
#include <modus/system/aspect.hpp>
#include <modus/kinematics/contact_constraint.hpp>

namespace modus
{

// This class computes kinematic contact constraints, a.k.a. normal and tangent
// velocity constraints.
class ContactKinematics : public Aspect {
 protected:
  System* system_;

 public:
  ContactKinematics(System* system);

  size_t NumNormalVelocityConstraints();

  size_t NumTangentVelocityConstraints();

  ContactConstraints GetNormalVelocityConstraints(const std::string& cs_mode, 
                                                  bool filter_contacts=false);

  ContactConstraints GetTangentVelocityConstraints(const std::string& cs_mode, 
                                                   const std::string& ss_mode,
                                                   bool filter_contacts=false);
};
MODUS_DEFINE_SHARED(ContactKinematics);

}