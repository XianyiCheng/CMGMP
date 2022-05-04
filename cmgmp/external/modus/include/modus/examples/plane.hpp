#pragma once
#include <modus/system/system.hpp>
#include <modus/system/body.hpp>
#include <modus/system/state.hpp>


namespace modus {
namespace examples {

Body* AddPlane(System* system, const Eigen::Vector3d& normal, double offset);

}
}