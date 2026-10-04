#include "task2_solver.hpp"

#include <print>
#include <string_view>
#include <vector>

void print_node_table(std::size_t n) {
    const double h   = 1.0 / static_cast<double>(n);
    const auto v     = solve_task2(n);
    const auto nodes = collect_node_data(v, h);

    std::println("\n--- Grid nodes check (n = {}) ---", n);
    std::println(
        "{:>4}  {:>8}  {:>16}  {:>16}  {:>14}  {:>14}",
        "i",
        "x_i",
        "v_calc",
        "u_exact",
        "|u - v|",
        "residual");

    for (const auto &node : nodes) {
        std::println(
            "{:>4}  {:>8.4f}  {:>16.10f}  {:>16.10f}  {:>14.6e}  {:>14.6e}",
            node.index,
            node.x,
            node.v_numeric,
            node.u_exact,
            node.error,
            node.residual);
    }
}

void print_experiment(const std::vector<std::size_t> &n_values, std::string_view label) {
    std::println("\n--- Experiment: {} ---", label);
    std::println(
        "{:>10}  {:>12}  {:>10}  {:>16}  {:>16}",
        "n",
        "step_h",
        "time_ms",
        "max_diff",
        "max_residual");

    for (std::size_t n : n_values) {
        const auto res = run_experiment(n);
        std::println(
            "{:>10}  {:>12.4e}  {:>10.4f}  {:>16.6e}  {:>16.6e}",
            res.n,
            res.h,
            res.elapsed_time_ms,
            res.max_error,
            res.max_residual);
    }
}

auto main() -> int {
    // 1. Detailed table for demonstration (n = 10)
    print_node_table(10);

    // 2. Experiment: Powers of 10
    const std::vector<std::size_t> powers_of_10 = {10, 100, 1'000, 10'000, 100'000, 1'000'000};
    print_experiment(powers_of_10, "Powers of 10");

    // 3. Experiment: Powers of 2
    const std::vector<std::size_t> powers_of_2 = {
        2,    4,    8,    16,    32,    64,    128,    256,    512,    1024,
        2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576};
    print_experiment(powers_of_2, "Powers of 2");

    return 0;
}
