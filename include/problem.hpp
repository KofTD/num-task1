#pragma once
#include <concepts>
#include <functional>

template <std::floating_point F>
struct Problem {
  std::function<F(F)> cstcc;  // cross section thermal conductivity coefficient
  std::function<F(F)> ihee;   // intensity of heat exchange environment
  std::function<F(F)> csds;   // cross‑section density of sinks

  // NOLINTBEGIN(readability-identifier-length)
  auto k(F x) -> F { return cstcc(x); }
  auto q(F x) -> F { return ihee(x); }
  auto f(F x) -> F { return csds(x); }
  // NOLINTEND(readability-identifier-length)
};
