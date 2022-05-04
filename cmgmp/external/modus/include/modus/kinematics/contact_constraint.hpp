#pragma once
#include <modus/common/eigen.hpp>
#include <modus/collision/contact.hpp>
#include <vector>
#include <memory>

namespace modus
{

struct ContactConstraint {
  enum Type {
    CC_NORMAL = 1,
    CC_TANGENT_X = 2,
    CC_TANGENT_Y = 3,
  };

  // A contact constraint is of the form lb <= nx <= ub.
  Eigen::RowVectorXd n;
  double lb;
  double ub;

  // Type of constraint.
  Type type;

  // Source contact which generated the constraint.
  ContactPtr source;
};

std::ostream& operator<<(std::ostream& os, const ContactConstraint& c);

using ContactConstraintPtr = std::shared_ptr<ContactConstraint>;
using ContactConstraints = std::vector<ContactConstraintPtr>;

std::ostream& operator<<(std::ostream& os, const ContactConstraints& C);

}