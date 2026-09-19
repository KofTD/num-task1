#pragma once
#include <algorithm>
#include <concepts>
#include <ranges>
#include <vector>

#include "fp-comprasion.hpp"
#include "grid.hpp"
#include "problem.hpp"

template <std::floating_point F>
struct TridiagonalMatrix {
  std::vector<F> a;
  std::vector<F> b;
  std::vector<F> c;
  F kappa1 = 0.0;
  F kappa2 = 0.0;

  TridiagonalMatrix() = default;
  TridiagonalMatrix(const TridiagonalMatrix&) = default;
  TridiagonalMatrix(TridiagonalMatrix&&) = default;
  auto operator=(const TridiagonalMatrix&) -> TridiagonalMatrix& = default;
  auto operator=(TridiagonalMatrix&&) -> TridiagonalMatrix& = default;
  explicit TridiagonalMatrix(std::integral auto size)
      : a(size - 2, 0.0), b(size - 2, 0.0), c(size - 2, 0.0) {
    if (size - 2 <= 0) {
      throw std::invalid_argument(
          "Size must be a positive number and greater than 2");
    }
  }

  [[nodiscard]] auto size() const -> std::size_t { return a.size() + 2; }
  ~TridiagonalMatrix() = default;
};

// Tridiagonal System of Algebraic Equations - TSAE
template <std::floating_point F>
struct TSAE {
  TridiagonalMatrix<F> matrix;
  std::vector<F> rhs;

  TSAE() = default;
  TSAE(const TSAE&) = default;
  TSAE(TSAE&&) = default;
  auto operator=(const TSAE&) -> TSAE& = default;
  auto operator=(TSAE&&) -> TSAE& = default;
  TSAE(Problem<F> problem, UniformGrid<F> grid) : TSAE(grid.nodes() + 1) {
    for (auto i : std::views::iota(static_cast<size_t>(1), grid.steps())) {
      const auto ai = problem.k(midpoint(grid.region(i - 1)));
      const auto ai1 = problem.k(midpoint(grid.region(i)));
      matrix.a[i - 1] = ai / (grid.step() * grid.step());
      matrix.b[i - 1] = ai1 / (grid.step() * grid.step());
      matrix.c[i - 1] =
          ((ai + ai1) / (grid.step() * grid.step())) + problem.q(grid.node(i));
      rhs[i] = problem.f(grid.node(i));
    }
  }
  explicit TSAE(std::unsigned_integral auto size)
      : matrix(size), rhs(size, 0.0) {}

  auto mu1(this auto&& self) -> auto&& { return self.rhs[0]; }
  auto mu2(this auto&& self) -> auto&& { return self.rhs.back(); }

  auto phi(this auto&& self, std::size_t index) -> auto&& {
    return self.rhs[index];
  }

  ~TSAE() = default;
};

template <std::floating_point F>
auto canUseTdma(const TSAE<F>& system) -> bool {
  const auto non_zero = [](const auto val) -> bool {
    return !approxEqual(std::abs(val), 0.0F);
  };
  // NOLINTBEGIN(readability-identifier-length)
  const auto non_strict_predominance = [](const auto& values) -> bool {
    const auto& [a, b, c] = values;
    return std::abs(c) >= std::abs(a) + std::abs(b);
  };
  const auto strict_predominance = [](const auto& values) -> bool {
    const auto& [a, b, c] = values;
    return std::abs(c) > std::abs(a) + std::abs(b);
  };
  // NOLINTEND(readability-identifier-length)
  auto non_zero_a = std::ranges::all_of(system.matrix.a, non_zero);
  auto non_zero_b = std::ranges::all_of(system.matrix.b, non_zero);
  auto non_strict_diagonal_predominance = std::ranges::all_of(
      std::views::zip(system.matrix.a, system.matrix.b, system.matrix.c),
      non_strict_predominance);
  auto strict_diagonal_predominance = std::ranges::all_of(
      std::views::zip(system.matrix.a, system.matrix.b, system.matrix.c),
      strict_predominance);

  auto condition1 = non_zero_a and non_zero_b and
                    non_strict_diagonal_predominance and
                    std::abs(system.matrix.kappa1) <= 1 and
                    std::abs(system.matrix.kappa2) < 1;
  auto condition2 = non_zero_a and non_zero_b and
                    strict_diagonal_predominance and
                    std::abs(system.matrix.kappa1) <= 1 and
                    std::abs(system.matrix.kappa2) <= 1;
  return condition1 or condition2;
}

// TDMA - Tridiagonal matrix algorithm
// This is progonka
template <std::floating_point F>
auto tdma(const TSAE<F>& system) -> std::vector<F> {
  const std::size_t matrix_size = system.matrix.size();
  const auto& matrix = system.matrix;
  std::vector<F> alpha(matrix_size - 1, 0.0);
  std::vector<F> beta(matrix_size - 1, 0.0);

  // Forward pass
  alpha[0] = system.matrix.kappa1;
  beta[0] = system.mu1();
  for (auto i : std::views::iota(static_cast<size_t>(1), alpha.size())) {
    F denom = matrix.c[i - 1] - (matrix.a[i - 1] * alpha[i - 1]);
    alpha[i] = matrix.b[i - 1] / denom;
    beta[i] = ((matrix.a[i - 1] * beta[i - 1]) + system.phi(i)) / denom;
  }
  // Backward pass
  std::vector<F> unknowns(matrix_size, 0.0);  // It's an y vector
  unknowns.back() = (system.mu2() + (system.matrix.kappa2 * beta.back())) /
                    (1.0 - (system.matrix.kappa2 * alpha.back()));
  for (auto i :
       std::views::iota(static_cast<std::size_t>(0), unknowns.size() - 1) |
           std::views::reverse) {
    unknowns[i] = (alpha[i] * unknowns[i + 1]) + beta[i];
  }
  return unknowns;
}
