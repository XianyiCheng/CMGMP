#include <modus/kinematics/contact_constraint.hpp>
#include <modus/common/linear_algebra.hpp>


using namespace modus;

std::ostream& modus::operator<<(std::ostream& os, const ContactConstraint& c)
{
  std::string type;
  switch (c.type) {
    case ContactConstraint::CC_NORMAL:
      type = "n";
      break;
    case ContactConstraint::CC_TANGENT_X:
      type = "x";
      break;
    case ContactConstraint::CC_TANGENT_Y:
      type = "y";
      break;
    default:
      MODUS_ASSERT(false);
  }
  os << "h: {" << std::setw(4) << c.lb << " <= " << PrettyMatrix(c.n,"",6) 
     << " <= " << c.ub << "}," << " type: " << type << ", src: " 
     << PrettyMatrix(c.source->point_A.transpose(),"",3);// << " "
    //  << PrettyMatrix(c.direction.transpose(),"",3);
  return os;
}

std::ostream& modus::operator<<(std::ostream& os, const ContactConstraints& C)
{
  for (size_t i = 0; i < C.size(); i++) {
    os << *C[i];
    if (i < C.size() - 1) {
      os << std::endl;
    }
  }
  return os;
}