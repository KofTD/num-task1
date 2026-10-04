#include "exp_utility.hpp"

#include <ranges>

using std::vector;

namespace rg = std::ranges;
namespace vw = std::views;

auto nodeDifference(const UniformGrid<double>& grid,
                    const std::vector<double>& tdma_result,
                    const std::function<double(double)>& target)
    -> vector<double> {
  auto true_answers = calcTrueAnswers(grid, target);
  return vw::zip_transform(std::minus<>(), true_answers, tdma_result) |
         rg::to<vector<double>>();
}

auto calcDiscrepancy(std::tuple<double, double, double> triplet, double a,
                     double b, double c, double phi) -> double {
  const auto [left, mid, right] = triplet;
  return (a * left) - (c * mid) + (b * right) + phi;
}

auto calcDiscrepancy(const TSAE<double>& system,
                     const vector<double>& tdma_result) -> vector<double> {
  const auto triplets = tdma_result | vw::adjacent<3>;
  const auto phis = vw::iota(1UZ, tdma_result.size() - 1) |
                    vw::transform([&system](std::size_t index) -> double {
                      return system.phi(index);
                    });
  const auto calc =
      static_cast<double (*)(std::tuple<double, double, double>, double, double,
                             double, double)>(calcDiscrepancy);
  const auto interior = vw::zip_transform(
      calc, triplets, system.matrix.a, system.matrix.b, system.matrix.c, phis);

  vector<double> discrepancy(tdma_result.size(), 0.0);
  discrepancy.front() = tdma_result.front() - system.mu1();
  rg::copy(interior, discrepancy.begin() + 1);
  discrepancy.back() = tdma_result.back() - system.mu2();
  return discrepancy;
}

auto calcDiscrepancy(const vector<double>& tdma_result,
                     const std::function<double(double)>& calc_a,
                     const std::function<double(double)>& calc_b,
                     const std::function<double(double)>& calc_c,
                     const std::function<double(double)>& calc_phi, double mu1,
                     double mu2) -> vector<double> {
  const auto triplets = tdma_result | vw::adjacent<3>;
  const auto phis =
      vw::iota(1UZ, tdma_result.size() - 1) | vw::transform(calc_phi);
  const auto calc =
      static_cast<double (*)(std::tuple<double, double, double>, double, double,
                             double, double)>(calcDiscrepancy);
  const auto a = vw::iota(1UZ, tdma_result.size() - 1) | vw::transform(calc_a);
  const auto b = vw::iota(1UZ, tdma_result.size() - 1) | vw::transform(calc_b);
  const auto c = vw::iota(1UZ, tdma_result.size() - 1) | vw::transform(calc_c);
  const auto interior = vw::zip_transform(calc, triplets, a, b, c, phis);

  vector<double> discrepancy(tdma_result.size(), 0.0);
  discrepancy.front() = tdma_result.front() - mu1;
  rg::copy(interior, discrepancy.begin() + 1);
  discrepancy.back() = tdma_result.back() - mu2;
  return discrepancy;
}
