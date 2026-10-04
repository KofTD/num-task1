#include "task2_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>

double exact_u(double x) { return (10.0 + 90.0 * x * x); }

double canonical_phi(double x) { return (450.0 * x * x - 2110.0); }

double compute_max_error(const std::vector<double>& v, double h) {
  std::size_t n = v.size() - 1;
  double max_err = 0.0;
  double x_i;

  for (size_t i = 0; i <= n; i++) {
    x_i = i * h;
    max_err = std::max(std::abs(v[i] - exact_u(x_i)), max_err);
  }

  return max_err;
}

double compute_max_residual(const std::vector<double>& v, double h) {
  std::size_t n = v.size() - 1;
  double max_res = std::max(std::abs(v[0] - 10.0), std::abs(v[n] - 100.0));

  double x_i;
  double lt,
      rt;  // lt - left part of the equation, rt - right part of the equation

  for (size_t i = 1; i < n; i++) {
    double x_i = i * h;
    double lt =
        12.0 * (v[i - 1] - 2.0 * v[i] + v[i + 1]) / (h * h) - 5.0 * v[i];
    double rt = 2110.0 - 450.0 * x_i * x_i;
    double current_res = std::abs(lt - rt);

    if (current_res > max_res) {
      max_res = current_res;
    }
  }

  return max_res;
}

std::vector<NodeData> collect_node_data(const std::vector<double>& v,
                                        double h) {
  std::size_t total_nodes = v.size();
  std::size_t n = total_nodes - 1;

  std::vector<NodeData> data;
  data.reserve(total_nodes);

  for (std::size_t i = 0; i <= n; ++i) {
    double x = i * h;
    double v_num = v[i];
    double u_ex = exact_u(x);
    double err = std::abs(u_ex - v_num);

    double res = 0.0;

    if (i == 0) {
      // v_0 = 10
      res = std::abs(v_num - 10.0);
    } else if (i == n) {
      // v_n = 100
      res = std::abs(v_num - 100.0);
    } else {
      double diff2 = (v[i - 1] - 2.0 * v[i] + v[i + 1]) / (h * h);
      double lt = 12.0 * diff2 - 5.0 * v[i];
      double rt = 2110.0 - 450.0 * x * x;
      res = std::abs(lt - rt);
    }

    data.push_back({.index = i,
                    .x = x,
                    .v_numeric = v_num,
                    .u_exact = u_ex,
                    .error = err,
                    .residual = res});
  }

  return data;
}

ExperimentResult run_experiment(std::size_t n) {
  auto start = std::chrono::high_resolution_clock::now();
  auto v = solve_task2(n);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsed = end - start;

  double h = 1.0 / n;
  double max_err = compute_max_error(v, h);
  double max_resid = compute_max_residual(v, h);

  ExperimentResult result = {n, h, elapsed.count(), max_err, max_resid};
  return result;
}
