#pragma once
#include <modus/common/memory.hpp>
#include <modus/system/aspect.hpp>


namespace modus
{

class CoefficientOfFriction : public Aspect {
 protected:
  double mu_;

 public:
  CoefficientOfFriction(double mu) { mu_ = mu; }

  void SetMu(double mu) { mu_ = mu; }

  double GetMu() { return mu_; }
};
using CoF = CoefficientOfFriction;
MODUS_DEFINE_SHARED(CoefficientOfFriction);
MODUS_DEFINE_SHARED(CoF);

}