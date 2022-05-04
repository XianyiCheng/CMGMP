#pragma once
#include <modus/system/composite.hpp>

namespace modus
{

template <typename T>
bool Composite::Has() const
{
  return (Get<T>() != nullptr);
}

template <typename T>
T* Composite::Get()
{
  AspectMap::iterator it = aspects_.find( typeid(T) );
  if (aspects_.end() == it) {
    return nullptr;
  }
  return static_cast<T*>(it->second.get());
}

// template <typename T>
// const T* Composite::Get() const
// {
//   return const_cast<Composite*>(this)->Get<T>();
// }

// template <typename B, typename T, typename... Args>
// T* Composite::Create(Args&&... args)
// {
//   T* aspect = new T(std::forward<Args>(args)...);
//   aspects_[typeid(B)] = std::unique_ptr<T>(aspect);

//   return aspect;
// }

template <typename B, typename T>
T* Composite::Create(T* aspect)
{
  // T* aspect = new T(std::forward<Args>(args)...);
  aspects_[typeid(B)] = std::unique_ptr<T>(aspect);

  return aspect;
}

}