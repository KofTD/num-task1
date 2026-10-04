#include "exp_utility.hpp"

#include <ranges>
#include <cstdlib>
#include <fstream>
#include <iostream>

using std::vector;

namespace rg = std::ranges;
namespace vw = std::views;

auto nodeDifference(
    const UniformGrid<double> &grid,
    const std::vector<double> &tdma_result,
    const std::function<double(double)> &target) -> vector<double> {
    auto true_answers = calcTrueAnswers(grid, target);
    return vw::zip_transform(std::minus<>(), true_answers, tdma_result) | rg::to<vector<double>>();
}

auto calcDiscrepancy(
    std::tuple<double, double, double> triplet, double a, double b, double c, double phi)
    -> double {
    const auto [left, mid, right] = triplet;
    return (a * left) - (c * mid) + (b * right) + phi;
}

auto calcDiscrepancy(const TSAE<double> &system, const vector<double> &tdma_result)
    -> vector<double> {
    const auto triplets = tdma_result | vw::adjacent<3>;
    const auto phis     = vw::iota(1UZ, tdma_result.size() - 1)
                        | vw::transform([&system](std::size_t index) -> double {
                          return system.phi(index);
                          });
    const auto calc =
        static_cast<double (*)(std::tuple<double, double, double>, double, double, double, double)>(
            calcDiscrepancy);
    const auto interior =
        vw::zip_transform(calc, triplets, system.matrix.a, system.matrix.b, system.matrix.c, phis);

    vector<double> discrepancy(tdma_result.size(), 0.0);
    discrepancy.front() = tdma_result.front() - system.mu1();
    rg::copy(interior, discrepancy.begin() + 1);
    discrepancy.back() = tdma_result.back() - system.mu2();
    return discrepancy;
}

auto calcDiscrepancy(
    const vector<double> &tdma_result,
    const UniformGrid<double> &grid,
    double a,
    double b,
    double c,
    const std::function<double(double)> &calc_phi,
    double mu1,
    double mu2) -> vector<double> {
    const auto triplets = tdma_result | vw::adjacent<3>;
    const auto nodes    = grid.nodeRange() | vw::drop(1) | vw::take(grid.nodes() - 2);
    const auto phis     = nodes | vw::transform(calc_phi);
    const auto calc     = [&](std::tuple<double, double, double> triplet, double phi) -> double {
        return calcDiscrepancy(triplet, a, b, c, phi);
    };
    const auto interior = vw::zip_transform(calc, triplets, phis);

    vector<double> discrepancy(tdma_result.size(), 0.0);
    discrepancy.front() = tdma_result.front() - mu1;
    rg::copy(interior, discrepancy.begin() + 1);
    discrepancy.back() = tdma_result.back() - mu2;
    return discrepancy;
}

auto saveResultsToFile(
    const UniformGrid<double>& grid, const Run& run, const std::string& filename) -> void {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Could not open file " << filename << " for writing\n";
        return;
    }
    out << "# x\tu(x)\tv\tdiff\tdiscrepancy\n";
    std::size_t i = 0;
    for (double x : grid.nodeRange()) {
        out << x << '\t' << run.exact[i] << '\t' << run.tdma[i] << '\t' << run.difference[i]
            << '\t' << run.discrepancy[i] << '\n';
        ++i;
    }
}

auto plotWithGnuplot(const std::string& data_filename, const std::string& title_prefix) -> void {
    const std::string script_filename = data_filename + ".gp";
    std::ofstream script(script_filename);
    if (!script.is_open()) {
        std::cerr << "Could not open gnuplot script file\n";
        return;
    }
    const std::string solution_png = data_filename + "_solution.png";
    const std::string error_png = data_filename + "_error.png";

    script << "set terminal pngcairo size 900,600\n";
    script << "set grid\n";
    script << "set output '" << solution_png << "'\n";
    script << "set title '" << title_prefix << ": u(x) vs v'\n";
    script << "set xlabel 'x'\n";
    script << "set ylabel 'u(x)'\n";
    script << "set key top left\n";
    script << "plot '" << data_filename << "' using 1:2 with lines lw 2 title 'u(x) exact', \\\n";
    script << "     '" << data_filename
        << "' using 1:3 with points pt 7 ps 0.5 title 'v (computed)'\n";
    script << "set output '" << error_png << "'\n";
    script << "set title '" << title_prefix << ": u(x) - v'\n";
    script << "set ylabel 'error'\n";
    script << "unset key\n";
    script << "plot '" << data_filename << "' using 1:4 with lines lw 2 title 'u(x)-v'\n";
    script.close();

    const std::string command = "gnuplot \"" + script_filename + "\"";
    int ret = std::system(command.c_str());
    if (ret != 0) {
        std::cerr << "gnuplot exited with code " << ret << " - is gnuplot installed?\n";
    }
}