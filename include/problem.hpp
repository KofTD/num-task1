#pragma once
#include <concepts>
#include <functional>

#include "grid.hpp"
#include "tsae.hpp"

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

template <std::floating_point F>
auto makeTSAE(Problem<F> problem, UniformGrid<F> grid) -> TSAE<F> {
  TSAE<F> tsae(grid.steps());
  auto& matrix = tsae.matrix;

  for (auto [i, region] : grid | std::views::enumerate) {
    auto ai = problem.k(midpoint(region));
    auto ai1 = problem.k(midpoint(grid.region(i)));
    matrix.a[i] = ai / (grid.step() * grid.step());
    matrix.b[i] = ai1 / (grid.step() * grid.step());
    matrix.c[i] = ((ai + ai1) / (grid.step() * grid.step())) +
                  problem.q(midpoint(region));
    tsae.rhs[i] = problem.f(midpoint(region));
  }
  return tsae;
}
