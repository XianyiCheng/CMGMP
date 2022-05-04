#pragma once
#include <memory>
#include <map>
#include <vector>
#include <modus/common/eigen.hpp>
#include <modus/system/body.hpp>
#include <modus/system/link.hpp>


namespace modus
{

struct Contact {
  // Public variables.
  // Body* body_A;
  // Body* body_B;
  // Link* link_A;
  // Link* link_B;
  int body_A;
  int body_B;
  int link_A;
  int link_B;
  Eigen::Vector3d point_A;
  Eigen::Vector3d point_B;
  Eigen::Vector3d contact_normal_B;
  double contact_distance;
  Eigen::Vector3d tangent_x;
  Eigen::Vector3d tangent_y;

  // For preprocessing use.
  int point_id_;
  int manifold_id_;
  double target_offset_;

  void SwapAB();
};

std::ostream& operator<<(std::ostream& os, const Contact& c);

using ContactPtr = std::shared_ptr<Contact>;
using Contacts = std::vector<ContactPtr>;
using ContactGroups = std::vector<Contacts>;
using ContactGroupMap = std::map<int, Contacts>;

}