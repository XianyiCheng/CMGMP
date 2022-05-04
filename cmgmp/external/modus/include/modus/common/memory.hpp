#pragma once
#include <memory>

#define MODUS_DEFINE_SHARED(X) using X ## Ptr = std::shared_ptr<X>