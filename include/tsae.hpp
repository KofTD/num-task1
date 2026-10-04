#pragma once
#include <algorithm>
#include <concepts>
#include <ranges>
#include <vector>

#include "fp-comparison.hpp"
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
  // auto row(std::size_t index) -> std::vector<F> {
  //   auto matrix_size = this->size();
  //   if (index < 0 or index > matrix_size - 1) {
  //     throw std::invalid_argument(
  //         std::format("Matrix has rows with numbers in [0, {}] not with {}",
  //                     matrix_size - 1, index));
  //   }

  //   std::vector<F> row(matrix_size, 0.0);
  //   if (index == 0) {
  //     row[0] = 1.0;
  //     row[1] = this->kappa1;
  //     return row;
  //   }
  //   if (index == this->size()) {
  //     row[index - 1] = 1.0;
  //     row[index - 2] = this->kappa2;
  //     return row;
  //   }
  //   std::size_t offset = index - 1;
  //   row[offset] = this->a[index];
  //   row[offset + 1] = this->c[index];
  //   row[offset + 2] = this->b[index];
  //   return row;
  // }

  // template <typename R>
  //   requires std::ranges::sized_range<R>
  // auto operator*(const R& vec) -> std::vector<F> {
  //   auto vec_size = std::size(vec);
  //   auto matrix_size = this->size();
  //   if (vec_size != matrix_size) {
  //     throw std::invalid_argument(
  //         std::format("Matrix and vector have incompatible sizes: {} vs {}",
  //                     matrix_size, vec_size));
  //   }
  //   std::vector<F> result;
  //   result.reserve(vec_size);
  //   auto rows =
  //       std::views::iota(0UZ, matrix_size) |
  //       std::views::transform([*this](std::size_t index) -> std::vector<F> {
  //         return this->row(index);
  //       });
  //   for (auto row : rows) {
  //     auto production = std::views::zip(row, vec) |
  //                       std::views::transform(std::multiplies<F>());
  //     auto res = std::ranges::fold_left(production, 0.0, std::plus<F>());
  //     result.push_back(res);
  //   }
  //   return result;
  // }

  [[nodiscard("result ignored")]] auto size() const -> std::size_t {
    return a.size() + 2;
  }
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
  TSAE(Problem<F> problem, UniformGrid<F> grid) : TSAE(grid.nodes()) {
    for (auto i : std::views::iota(1UZ, grid.steps())) {
      // NOLINTBEGIN(*-identifier-length)
      const auto a = problem.k(midpoint(grid.region(i - 1)));
      const auto a_next = problem.k(midpoint(grid.region(i)));
      const auto d = problem.q(grid.node(i));
      const auto phi = problem.f(grid.node(i));
      const auto step_sqr = grid.step() * grid.step();
      matrix.a[i - 1] = a / step_sqr;
      matrix.b[i - 1] = a_next / step_sqr;
      matrix.c[i - 1] = ((a + a_next) / step_sqr) + d;
      rhs[i] = phi;
      // NOLINTEND(*-identifier-length)
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
    return !approxEqual(std::abs(val), 0.0);
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
  for (auto i : std::views::iota(1UZ, alpha.size())) {
    F denom = matrix.c[i - 1] - (matrix.a[i - 1] * alpha[i - 1]);
    alpha[i] = matrix.b[i - 1] / denom;
    beta[i] = ((matrix.a[i - 1] * beta[i - 1]) + system.phi(i)) / denom;
  }
  // Backward pass
  std::vector<F> unknowns(matrix_size, 0.0);  // It's an y vector
  unknowns.back() = (system.mu2() + (system.matrix.kappa2 * beta.back())) /
                    (1.0 - (system.matrix.kappa2 * alpha.back()));
  for (auto i :
       std::views::iota(0UZ, unknowns.size() - 1) | std::views::reverse) {
    unknowns[i] = (alpha[i] * unknowns[i + 1]) + beta[i];
  }
  return unknowns;
}
