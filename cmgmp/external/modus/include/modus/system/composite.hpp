#pragma once
#include <map>
#include <memory>
#include <typeinfo>
#include <typeindex>
#include <modus/common/memory.hpp>
#include <modus/system/aspect.hpp>


namespace modus
{

class Composite : public std::enable_shared_from_this<Composite> {
 protected:
  using AspectMap = std::map<std::type_index, std::shared_ptr<Aspect> >;
  AspectMap aspects_;

 public:
  virtual ~Composite() = default;
  Composite() = default;

  Composite(const Composite&) = delete;

  Composite(Composite&&) = delete;

  Composite& operator=(const Composite&) = delete;

  Composite& operator=(Composite&&) = delete;

  template <typename T>
  bool Has() const;

  template <typename T>
  T* Get();

  // template <typename T>
  // const T* Get() const;

  // template <typename B, typename T, typename... Args>
  // T* Create(Args&&... args);

  template <typename B, typename T>
  T* Create(T* aspect);
};
MODUS_DEFINE_SHARED(Composite);

}

#include <modus/system/detail/composite.hpp>