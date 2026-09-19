#pragma once
#include <algorithm>
#include <concepts>
#include <ranges>
#include <vector>

template <std::floating_point F>
struct TridiagonalMatrix {
  std::vector<F> a;
  std::vector<F> b;
  std::vector<F> c;
  F kappa1 = 0.0;
  F kappa2 = 0.0;

  TridiagonalMatrix() = default;
  explicit TridiagonalMatrix(std::size_t size)
      : a(size - 2, 0.0), b(size - 2, 0.0), c(size - 2, 0.0) {}

  [[nodiscard]] auto size() const -> std::size_t { return a.size() + 2; }
};

// Tridiagonal System of Algebraic Equations - TSAE
template <std::floating_point F>
struct TSAE {
  TridiagonalMatrix<F> matrix;
  std::vector<F> rhs;

  auto mu1(this auto&& self) -> auto&& { return self.rhs[0]; }
  auto mu2(this auto&& self) -> auto&& { return self.rhs.back(); }

  auto phi(this auto&& self, std::size_t index) -> auto&& {
    return self.rhs[index];
  }
};

template <std::floating_point F>
auto canUseTdma(const TSAE<F>& system) -> bool {
  const auto non_zero = [](const auto val) -> bool {
    return std::abs(val) != 0;
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
  const int32_t matrix_size = system.matrix.size();
  std::vector<F> alpha(matrix_size, 0.0);
  std::vector<F> beta(matrix_size, 0.0);
  const auto& matrix = system.matrix;

  // Forward pass
  alpha[0] = system.matrix.kappa1;
  beta[0] = system.mu1();
  for (auto i : std::views::iota(0, matrix_size - 1)) {
    F denom = matrix.c[i] - (matrix.a[i] * alpha[i]);
    alpha[i + 1] = matrix.b[i] / denom;
    beta[i + 1] = ((matrix.a[i] * beta[i]) + system.phi(i)) / denom;
  }

  // Backward pass
  std::vector<F> unknowns(matrix_size, 0.0);  // It's an y vector
  unknowns[matrix_size] =
      (system.mu2() + (system.matrix.kappa2 * beta[matrix_size])) /
      (1.0 - (system.matrix.kappa2 * alpha[matrix_size]));
  for (auto i : std::views::iota(0, matrix_size) | std::views::reverse) {
    unknowns[i] = (alpha[i + 1] * unknowns[i + 1]) + beta[i + 1];
  }
  return unknowns;
}
