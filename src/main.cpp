#include <cstdlib>
#include <print>
#include <ranges>
#include <vector>

#include "problem.hpp"
#include "tsae.hpp"

using std::vector;

// NOLINTBEGIN
template <std::floating_point F>
auto f(F x) -> F {
  return 450 * x * x - 2110;
}
template <std::floating_point F>
auto k(F x) -> F {
  return 12.0;
}
template <std::floating_point F>
auto q(F x) -> F {
  return 5.0;
}
template <std::floating_point F>
auto u(F x) -> F {
  return 10 + 90 * x * x;
}
// NOLINTEND

auto main() -> int {
  using fp_type = double;
  const Problem<fp_type> problem = {
      .cstcc = k<fp_type>,
      .ihee = q<fp_type>,
      .csds = f<fp_type>,
  };
  const UniformGrid<fp_type> grid(0.0, 1.0, 8);
  TSAE<fp_type> system(problem, grid);
  // NOLINTNEXTLINE(*-magic-numbers)
  system.mu1() = 10.0;
  system.mu2() = 100.0;
  system.matrix.kappa1 = 0.0;
  system.matrix.kappa2 = 0.0;
  const auto nodes =
      std::views::iota(0UZ, grid.nodes()) |
      std::views::transform(
          [&grid](std::size_t index) -> fp_type { return grid.node(index); });
  const auto standards = nodes | std::views::transform(u<fp_type>);
  const auto results = tdma(system);
  for (auto [point, res, standard] :
       std::views::zip(nodes, results, standards)) {
    std::println("u'({0}) = {1}, u({0}) = {2} ", point, res, standard);
  }
  // } else {
  //   std::println("Can not use TDMA");
  // }

  return 0;
}
