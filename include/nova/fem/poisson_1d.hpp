#pragma once

#include "nova/linalg/conjugate_gradient.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

struct Poisson1DProblem {
    double length{1.0};
    std::size_t elements{20};
    double conductivity{1.0};
    std::function<double(double)> source;
    double left_value{0.0};
    double right_value{0.0};
};

struct Poisson1DResult {
    std::vector<double> coordinates;
    std::vector<double> values;
    linalg::SolverResult linear_solver;
};

// Solves -(k*u')' = f on [0, length] with prescribed endpoint values.
// Uses uniform linear P1 elements and two-point Gauss quadrature.
[[nodiscard]] Poisson1DResult solve_poisson_1d(
    const Poisson1DProblem& problem,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::fem
