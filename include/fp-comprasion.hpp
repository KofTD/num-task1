#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdlib>
#include <limits>

template <std::floating_point F>
constexpr auto approxEqual(
    F a, F b, F relative_tolerance = std::numeric_limits<F>::epsilon() * 7,
    F absolute_tolerance = std::numeric_limits<F>::epsilon()) -> bool {
  if (a == b) {
    return true;
  }
  const F diff = std::abs(a - b);
  return diff <= absolute_tolerance or
         diff <= relative_tolerance * std::max(std::abs(a), std::abs(b));
}
