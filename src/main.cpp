#include "exp_utility.hpp"
#include "optimized-tdma.hpp"
#include "problem.hpp"
#include "tsae.hpp"

#include <boost/program_options.hpp>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <print>
#include <ranges>
#include <vector>

using std::vector;
using u32 = std::uint32_t;
using std::chrono::steady_clock;

namespace po     = boost::program_options;
namespace views  = std::views;
namespace ranges = std::ranges;

// NOLINTBEGIN
template <std::floating_point F>
auto f(F x) -> F {
    return 450 * x * x - 2110;
}

template <std::floating_point F>
auto k(F x) -> F {
    return 12.0;
}

template <std::floating_point F>
auto q(F x) -> F {
    return 5.0;
}

template <std::floating_point F>
auto u(F x) -> F {
    return 10 + 90 * x * x;
}

// NOLINTEND

enum class RunType : std::int8_t {
    Optimized,
    Universal,
};

auto universalRun(const Problem<double> &problem, const UniformGrid<double> &grid) -> Run {
    auto start = steady_clock::now();
    TSAE<double> system(problem, grid);
    system.mu1()         = 10.0;
    system.mu2()         = 100.0;
    system.matrix.kappa1 = 0.0;
    system.matrix.kappa2 = 0.0;
    auto result          = tdma(system);
    auto finish          = steady_clock::now();
    auto time            = finish - start;
    auto node_difference = nodeDifference(grid, result, u<double>);
    auto discrepancy     = calcDiscrepancy(system, result);
    auto true_answers    = calcTrueAnswers(grid, u<double>) | std::ranges::to<vector>();
    return {
        .exact       = std::move(true_answers),
        .tdma        = std::move(result),
        .difference  = std::move(node_difference),
        .discrepancy = std::move(discrepancy),
        .run_time    = time,
    };
}

auto prepareForRun(u32 number_of_nodes) -> std::pair<Problem<double>, UniformGrid<double>> {
    Problem<double> problem = {
        .cstcc = k<double>,
        .ihee  = q<double>,
        .csds  = f<double>,
    };
    UniformGrid<double> grid(0.0, 1.0, number_of_nodes);
    return {problem, grid};
}

auto optimizedRun(const UniformGrid<double> &grid) -> Run {
    auto start           = steady_clock::now();
    auto result          = optimizedTdma(grid);
    auto finish          = steady_clock::now();
    auto time            = finish - start;
    auto node_difference = nodeDifference(grid, result, u<double>);
    auto discrepancy     = calcDiscrepancy(
        result,
        grid,
        calcA(grid.step()),
        calcB(grid.step()),
        calcC(grid.step()),
        calcPhi,
        10.0,
        100.0);
    auto true_answers = calcTrueAnswers(grid, u<double>) | std::ranges::to<vector>();
    return {
        .exact       = std::move(true_answers),
        .tdma        = std::move(result),
        .difference  = std::move(node_difference),
        .discrepancy = std::move(discrepancy),
        .run_time    = time,
    };
}

auto printTable1(RunType rtype, std::size_t number_of_nodes) -> void {
    auto [problem, grid] = prepareForRun(number_of_nodes);
    Run res{};
    std::string srtype;
    switch (rtype) {
        case RunType::Universal:
            res    = universalRun(problem, grid);
            srtype = "universal";
            break;
        case RunType::Optimized:
            res    = optimizedRun(grid);
            srtype = "optimized";
            break;
    }
    std::println("Table 1. Implementation: {}, n = {}", srtype, number_of_nodes);
    std::println(
        "{:<14} {:^14} {:^20} {:^20} {:^16} {:^14}",
        "Node number",
        "x",
        "u(x)",
        "v",
        "u(x) - v",
        "discrepancy");
    const auto row_data =
        std::views::zip(grid.nodeRange(), res.exact, res.tdma, res.difference, res.discrepancy);
    for (auto [num, row] : std::views::enumerate(row_data)) {
        auto [node, exact, tdma, difference, discrepancy] = row;
        std::println(
            "{:<14} {:^14.8f} {:^20.12f} {:^20.12f} {:^16.3e} {:^14.3e}",
            num,
            node,
            exact,
            tdma,
            difference,
            discrepancy);
    }
    const auto max_error       = std::ranges::max_element(res.difference, {}, std::abs);
    const auto max_discrepancy = std::ranges::max_element(res.discrepancy, {}, std::abs);
    std::println(
        "max |u(x) - v| = {:.6e} at {}",
        std::abs(*max_error),
        std::distance(res.difference.begin(), max_error));
    std::println(
        "max |discrepancy| = {:.6e} at {}",
        std::abs(*max_discrepancy),
        std::distance(res.discrepancy.begin(), max_discrepancy));
    std::println("time = {:.6e} s", std::chrono::duration<double>(res.run_time).count());
}

template <std::integral I1, std::integral I2, std::integral I = std::common_type_t<I1, I2>>
auto pow(I1 base, I2 pow) -> I {
    I ibase = base;
    I ipow  = pow;

    I result = 1;
    for ([[maybe_unused]] auto blk : std::views::iota(static_cast<I>(0), ipow)) {
        result *= ibase;
    }
    return result;
}

template <std::ranges::input_range R>
    requires std::convertible_to<std::ranges::range_value_t<R>, u32>
auto printTable2Rows(R &&numbers_of_nodes, RunType rtype) -> void {
    const auto max_by_abs_val = std::bind_back(
        std::ranges::max_element, std::less<>(), static_cast<double (*)(double)>(std::abs));
    std::string srtype;
    switch (rtype) {
        case RunType::Universal:
            srtype = "universal";
            break;
        case RunType::Optimized:
            srtype = "optimized";
            break;
    }
    std::println("Table 2: implementation {}", srtype);
    std::println(
        "{:<20} {:^16} {:^18} {:^12}",
        "Number of regions",
        "max error",
        "max discrepancy",
        "time, s");
    for (auto &&num : numbers_of_nodes) {
        const auto [problem, grid] = prepareForRun(num);
        Run run{};
        switch (rtype) {
            case RunType::Universal:
                run = universalRun(problem, grid);
                break;
            case RunType::Optimized:
                run = optimizedRun(grid);
                break;
        }
        const auto max_error       = *max_by_abs_val(run.difference);
        const auto max_discrepancy = *max_by_abs_val(run.discrepancy);
        const auto time            = std::chrono::duration<double>(run.run_time).count();
        std::println(
            "{:<20} {:^16.6e} {:^18.6e} {:^12.6e}",
            num,
            std::abs(max_error),
            std::abs(max_discrepancy),
            time);
    }
}

auto printTable2(RunType rtype) -> void {
    const auto pow10       = std::bind_front(pow<u32, u32>, 10);
    const auto pow2        = std::bind_front(pow<u32, u32>, 2);
    const auto upper_bound = 1'000'000;
    const u32 log2_up      = static_cast<u32>(std::ceil(std::log2(upper_bound) + 1));
    const u32 log10_up     = static_cast<u32>(std::ceil(std::log10(upper_bound) + 1));
    auto numbers_of_nodes  = std::views::iota(1UL, log10_up) | views::transform(pow10);
    printTable2Rows(numbers_of_nodes, rtype);
    numbers_of_nodes = std::views::iota(1UL, log2_up) | views::transform(pow2);
    printTable2Rows(numbers_of_nodes, rtype);
}

auto main(int argc, char *argv[]) -> int {
    std::string algo;
    u32 number_of_nodes = 0;
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "this help message")(
        "algorithm,a",
        po::value<std::string>(&algo)->required(),
        "choose version of algorithm: (universal|optimized)")(
        "nodes,n", po::value<u32>(&number_of_nodes)->required(), "number of nodes in grid");
    po::variables_map vmap;
    po::store(po::parse_command_line(argc, argv, desc), vmap);
    po::notify(vmap);

    if (vmap.contains("help") or vmap.empty()) {
        std::cout << desc << '\n';
        return 0;
    }

    RunType rtype = RunType::Universal;
    if (algo == "universal") {
        rtype = RunType::Universal;
    } else if (algo == "optimized") {
        rtype = RunType::Optimized;
    } else {
        std::println("Algorithm must be universal or optimized, not {}", algo);
        return 1;
    }

    printTable1(rtype, number_of_nodes);
    std::println("======================");
    printTable2(rtype);

    return 0;
}
