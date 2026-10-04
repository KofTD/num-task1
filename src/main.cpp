#include <boost/program_options.hpp>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <print>
#include <ranges>
#include <vector>

#include "exp_utility.hpp"
#include "problem.hpp"
#include "tsae.hpp"

using std::vector;
using u32 = std::uint32_t;
using std::chrono::steady_clock;

namespace po = boost::program_options;
namespace views = std::views;
namespace ranges = std::ranges;

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

auto universalRun(u32 number_of_nodes) -> Run {
  Problem<double> problem = {
      .cstcc = k<double>,
      .ihee = q<double>,
      .csds = f<double>,
  };
  UniformGrid<double> grid(0.0, 1.0, number_of_nodes);
  auto start = steady_clock::now();
  TSAE<double> system(problem, grid);
  system.mu1() = 10.0;
  system.mu2() = 100.0;
  system.matrix.kappa1 = 0.0;
  system.matrix.kappa2 = 0.0;
  auto result = tdma(system);
  auto finish = steady_clock::now();
  auto time = finish - start;
  auto node_difference = nodeDifference(grid, result, u<double>);
  auto discrepancy = calcDiscrepancy(system, result);
  auto true_answers =
      calcTrueAnswers(grid, u<double>) | std::ranges::to<vector>();
  return {
      .exact = std::move(true_answers),
      .tdma = std::move(result),
      .difference = std::move(node_difference),
      .discrepancy = std::move(discrepancy),
      .run_time = time,
  };
}

auto optimizedRun() -> void {}

auto printTable1(std::size_t number_of_nodes) {
  const Run res = universalRun(number_of_nodes);
  std::println("Table 1. Implementation: {}, n = {}", "universal",
               number_of_nodes);
  std::println("{:>8} {:>14} {:>20} {:>20} {:>12} {:>12}", "i", "x_i",
               "u_i (exact)", "v_i (tdma)", "u_i - v_i", "r_i");
  const double h = 1.0 / number_of_nodes;
  for (std::size_t i = 0; i <= number_of_nodes; ++i) {
    std::println("{:>8} {:>14.8f} {:>20.12f} {:>20.12f} {:>12.3e} {:>12.3e}", i,
                 static_cast<double>(i) * h, res.exact[i], res.tdma[i],
                 res.difference[i], res.discrepancy[i]);
  }
  std::println(
      "max |u_i - v_i| = {:.6e}",
      *std::ranges::max_element(res.difference, std::less<>(),
                                static_cast<double (*)(double)>(std::abs)));
  std::println(
      "max |r_i|       = {:.6e}",
      *std::ranges::max_element(res.discrepancy, std::less<>(),
                                static_cast<double (*)(double)>(std::abs)));
  std::println("time            = {:.6e} s",
               std::chrono::duration<double>(res.run_time).count());
}
auto main(int argc, char* argv[]) -> int {
  std::string algo;
  u32 number_of_nodes = 0;
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "this help message")(
      "algorithm,a", po::value<std::string>(&algo)->required(),
      "choose version of algorithm: (universal|optimized)")(
      "nodes,n", po::value<u32>(&number_of_nodes)->required(),
      "number of nodes in grid");
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);

  if (vmap.contains("help") or vmap.empty()) {
    std::cout << desc << '\n';
    return 0;
  }

  if (algo != "universal" and algo != "optimized") {
    std::println("Algorithm must be universal or optimized, not {}", algo);
    return 1;
  }

  if (algo == "universal") {
    printTable1(number_of_nodes);
  } else {
    optimizedRun();
  }

  return 0;
}
