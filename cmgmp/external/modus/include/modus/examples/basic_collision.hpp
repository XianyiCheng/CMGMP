#pragma once
#include <modus/system/system.hpp>
#include <modus/collision/collision.hpp>
#include <modus/geometry/box.hpp>
#include <modus/geometry/plane.hpp>

namespace modus {
namespace examples {

class BasicCollision : public Collision {
 protected:
  System* system_;
  double  margin_;

 public:
  BasicCollision(System* system);

  void SetMargin(double margin) { margin_ = margin; }

  Contacts GetContactPoints() override;

  Contacts Collide(BoxGeometry* box, PlaneGeometry* plane);
};

}
}