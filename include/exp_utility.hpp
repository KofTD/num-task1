#pragma once
#include "grid.hpp"
#include "tsae.hpp"

#include <chrono>
#include <vector>

struct Run {
    std::vector<double> exact;
    std::vector<double> tdma;
    std::vector<double> difference;
    std::vector<double> discrepancy;
    std::chrono::nanoseconds run_time;
};

inline auto
calcTrueAnswers(const UniformGrid<double> &grid, const std::function<double(double)> &target) {
    auto nodes  = grid.nodeRange();
    auto values = std::views::transform(nodes, target);
    return values;
}

auto nodeDifference(
    const UniformGrid<double> &grid,
    const std::vector<double> &tdma_result,
    const std::function<double(double)> &target) -> std::vector<double>;

auto calcDiscrepancy(const TSAE<double> &system, const std::vector<double> &tdma_result)
    -> std::vector<double>;
auto calcDiscrepancy(
    const std::vector<double> &tdma_result,
    const std::function<double(double)> &calc_a,
    const std::function<double(double)> &calc_b,
    const std::function<double(double)> &calc_c,
    const std::function<double(double)> &calc_phi,
    double mu1,
    double mu2) -> std::vector<double>;
