#pragma once

#include "grid.hpp"

#include <vector>

// NOLINTBEGIN(*-magic-numbers)
auto calcPhi(double x) -> double {
    return (450.0 * x * x) - 2110.0;
}

auto calcA(double step) -> double {
    return 12.0 / (step * step);
};

auto calcB(double step) -> double {
    return 12.0 / (step * step);
};

auto calcC(double step) -> double {
    return (24.0 / (step * step)) + 5.0;
};

// NOLINTEND(*-magic-numbers)

auto optimizedTdma(const UniformGrid<double> &grid) -> std::vector<double> {
    const auto step  = grid.step();
    const auto A     = calcA(step);
    const auto B     = calcB(step);
    const auto C     = calcC(step);
    const auto nodes = grid.nodes();

    std::vector<double> alpha;
    std::vector<double> beta;
    alpha.reserve(nodes - 1);
    beta.reserve(nodes - 1);

    // Forward pass
    alpha.push_back(0.0);
    beta.push_back(0.0);
    const auto phis = std::views::iota(1UZ, alpha.size())
                    | std::views::transform([step](std::size_t i) -> double {
                          return calcPhi(static_cast<double>(i) * step);
                      });
    for (const auto phi : phis) {
        const auto denom      = C - (A * alpha.back());
        const auto next_alpha = B / denom;
        const auto next_beta  = ((A * beta.back()) + phi) / denom;
        alpha.push_back(next_alpha);
        beta.push_back(next_beta);
    }

    // Backward pass
    std::vector<double> v(nodes, double{0});
    v.back()    = 100.0;
    double next = v.back();
    for (auto &&[al, be, out] :
         std::views::zip(alpha, beta, v | std::views::take(nodes - 1)) | std::views::reverse) {
        next = (al * next) + be;
        out  = next;
    }
    return v;
}
