#pragma once
#include <chrono>
#include <vector>

#include "grid.hpp"
#include "tsae.hpp"

struct Run {
  std::vector<double> exact;
  std::vector<double> tdma;
  std::vector<double> difference;
  std::vector<double> discrepancy;
  std::chrono::nanoseconds run_time;
};

inline auto calcTrueAnswers(const UniformGrid<double>& grid,
                            const std::function<double(double)>& target) {
  auto nodes = std::views::iota(0UZ, grid.nodes()) |
               std::views::transform([grid](std::size_t index) -> double {
                 return grid.node(index);
               });
  auto values = std::views::transform(nodes, target);
  return values;
}

auto nodeDifference(const UniformGrid<double>& grid,
                    const std::vector<double>& tdma_result,
                    const std::function<double(double)>& target)
    -> std::vector<double>;

auto calcDiscrepancy(const TSAE<double>& system,
                     const std::vector<double>& tdma_result)
    -> std::vector<double>;
