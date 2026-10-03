#include "exp_utility.hpp"

#include <ranges>

using std::vector;

namespace views = std::views;
namespace ranges = std::ranges;

auto nodeDifference(const UniformGrid<double>& grid,
                    const std::vector<double>& tdma_result,
                    const std::function<double(double)>& target)
    -> vector<double> {
  auto true_answers = calcTrueAnswers(grid, target);
  return views::zip_transform(std::minus<>(), true_answers, tdma_result) |
         ranges::to<vector<double>>();
}
auto calcDiscrepancy(const TSAE<double>& system,
                     const vector<double>& tdma_result) -> vector<double> {
  auto triplets = tdma_result | views::adjacent<3>;
  auto phis = views::iota(1UZ, tdma_result.size() - 1) |
              views::transform([&system](std::size_t index) -> double {
                return system.phi(index);
              });
  auto calc_r = [](auto triplet, double a, double b, double c, double phi) {
    auto [left, mid, right] = triplet;
    return (a * left) - (c * mid) + (b * right) + phi;
  };
  auto interior =
      std::views::zip_transform(calc_r, triplets, system.matrix.a,
                                system.matrix.b, system.matrix.c, phis);

  vector<double> discrepancy(tdma_result.size(), 0.0);
  discrepancy.front() = tdma_result.front() - system.mu1();
  std::ranges::copy(interior, discrepancy.begin() + 1);
  discrepancy.back() = tdma_result.back() - system.mu2();
  return discrepancy;
}
