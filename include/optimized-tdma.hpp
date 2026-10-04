#pragma once

#include <vector>

#include "grid.hpp"

// NOLINTBEGIN(*-magic-numbers)
auto calcPhi(double x) -> double { return (450.0 * x * x) - 2110.0; }
auto calcA(double step) -> double { return 12.0 / (step * step); };
auto calcB(double step) -> double { return 12.0 / (step * step); };
auto calcC(double step) -> double { return (24.0 / (step * step)) + 5.0; };
// NOLINTEND(*-magic-numbers)

auto optimizedTdma(const UniformGrid<double>& grid) -> std::vector<double> {
  auto step = grid.step();
  double A = calcA(step);
  double B = calcB(step);
  double C = calcC(step);

  std::vector<double> alpha(grid.nodes() - 1);
  std::vector<double> beta(grid.nodes() - 1);
  alpha.front() = 0.0;
  beta.front() = 10.0;
  double x_i = 0.0;
  double f_i = 0.0;
  double denom = 0.0;
  for (size_t i = 1; i < alpha.size(); i++) {
    x_i = i * step;
    f_i = calcPhi(x_i);
    denom = C - (A * alpha[i - 1]);

    alpha[i] = B / denom;
    beta[i] = ((A * beta[i - 1]) + f_i) / denom;
  }

  std::vector<double> v(grid.nodes());
  v.back() = 100.0;
  for (int i = v.size() - 2; i >= 0; i--) {
    v[i] = (alpha[i] * v[i + 1]) + beta[i];
  }
  return v;
}
