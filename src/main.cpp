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
#include "optimized-tdma.hpp"
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

auto universalRun(const Problem<double>& problem,
                  const UniformGrid<double>& grid) -> Run {
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

auto prepareForRun(u32 number_of_nodes)
    -> std::pair<Problem<double>, UniformGrid<double>> {
  Problem<double> problem = {
      .cstcc = k<double>,
      .ihee = q<double>,
      .csds = f<double>,
  };
  UniformGrid<double> grid(0.0, 1.0, number_of_nodes);
  return {problem, grid};
}

auto optimizedRun(const UniformGrid<double>& grid) -> Run {
  auto start = steady_clock::now();
  auto result = optimizedTdma(grid);
  auto finish = steady_clock::now();
  auto time = finish - start;
  auto node_difference = nodeDifference(grid, result, u<double>);
  auto discrepancy =
      calcDiscrepancy(result, calcA, calcB, calcC, calcPhi, 10.0, 100.0);
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

auto printTable1(std::string_view run_type, std::size_t number_of_nodes)
    -> void {
  auto [problem, grid] = prepareForRun(number_of_nodes);
  const Run res = run_type == "optimized" ? optimizedRun(grid)
                                          : universalRun(problem, grid);
  std::println("Table 1. Implementation: {}, n = {}",
               run_type == "optimized" ? "optimized" : "universal",
               number_of_nodes);
  std::println("{:>14} {:>14} {:>20} {:>20} {:>12} {:>14}", "Node number", "x",
               "u(x)", "v", "u(x) - v", "discrepancy");
  const auto nodes = UniformGrid<double>(0.0, 1.0, number_of_nodes).nodeRange();
  const auto row_data = std::views::zip(nodes, res.exact, res.tdma,
                                        res.difference, res.discrepancy);
  for (auto [num, row] : std::views::enumerate(row_data)) {
    auto [node, exact, tdma, difference, discrepancy] = row;
    std::println("{:>14} {:>14.8f} {:>20.12f} {:>20.12f} {:>12.3e} {:>14.3e}",
                 num, node, exact, tdma, difference, discrepancy);
  }
  const auto max_by_abs_val =
      std::bind_back(std::ranges::max_element, std::less<>(),
                     static_cast<double (*)(double)>(std::abs));
  const auto max_error = *max_by_abs_val(res.difference);
  const auto max_discrepancy = *max_by_abs_val(res.discrepancy);
  std::println("max |u(x) - v| = {:.6e}", std::abs(max_error));
  std::println("max |discrepancy| = {:.6e}", std::abs(max_discrepancy));
  std::println("time = {:.6e} s",
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

  printTable1(algo, number_of_nodes);

  return 0;
}
