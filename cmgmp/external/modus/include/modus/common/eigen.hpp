#pragma once
#include <modus/common/assert.hpp>
#define eigen_assert(x) \
  if (!(x)) { MODUS_ASSERT((x)); }
#include <Eigen/Dense>