#pragma once

#include <cstddef>
#include <vector>

// struct for table
struct NodeData {
  std::size_t index;  // index i
  double x;           // coord x_i
  double v_numeric;   // solved value v_i
  double u_exact;     // exact value u(x_i)
  double error;       // |u_exact - v_numeric|
  double residual;    // nevuazka r_i
};

// struct for the results of one computational experiment
struct ExperimentResult {
  std::size_t n;           // The number of grid sections
  double h;                // step (1.0 / n)
  double elapsed_time_ms;  // solving time
  double max_error;
  double max_residual;
};

// exact solve u(x) = 10 + 90 * x^2
double exact_u(double x);

// phi(x) = 450 * x^2 - 2110
double canonical_phi(double x);

// returns vector v_i, length (n + 1)
std::vector<double> solve_task2(std::size_t n);

double compute_max_error(const std::vector<double>& v, double h);

double compute_max_residual(const std::vector<double>& v, double h);

// collection of complete information for all nodes
std::vector<NodeData> collect_node_data(const std::vector<double>& v, double h);

ExperimentResult run_experiment(std::size_t n);