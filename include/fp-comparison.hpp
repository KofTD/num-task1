#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstdlib>
#include <limits>

template <std::floating_point F1, std::floating_point F2,
          std::floating_point F = std::common_type_t<F1, F2>>
constexpr auto approxEqual(
    // NOLINTNEXTLINE(*-identifier-length, *-magic-numbers)
    F1 a, F2 b, F relative_tolerance = std::numeric_limits<F>::epsilon() * 7,
    F absolute_tolerance = std::numeric_limits<F>::epsilon()) -> bool {
  F coma = a;
  F comb = b;
  if (coma == comb) {
    return true;
  }
  const F diff = std::abs(coma - comb);
  return diff <= absolute_tolerance or
         diff <= relative_tolerance * std::max(std::abs(coma), std::abs(comb));
}
