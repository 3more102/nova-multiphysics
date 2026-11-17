#pragma once

#include "nova/nonlinear/newton_solver.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace nova::fem {

struct NonlinearDiffusion1DProblem {
    double length{1.0};
    std::size_t elements{32};
    std::function<double(double)> source;
    double left_value{0.0};
    double right_value{0.0};
};

struct NonlinearDiffusion1DResult {
    std::vector<double> coordinates;
    nonlinear::NewtonResult nonlinear_solver;
};

// Solves -(d/dx[(1+u^2) du/dx]) = f(x) on a uniform P1 mesh.
// Uses three-point Gauss quadrature and strong endpoint Dirichlet conditions.
// An empty initial guess selects zero at all nodes; endpoint values are imposed.
[[nodiscard]] NonlinearDiffusion1DResult solve_nonlinear_diffusion_1d(
    const NonlinearDiffusion1DProblem& problem,
    std::span<const double> initial_guess = {},
    const nonlinear::NewtonOptions& newton_options = {}
);

} // namespace nova::fem
