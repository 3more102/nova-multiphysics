#pragma once

#include "nova/fem/h1_space_2d.hpp"
#include "nova/fem/mesh2d.hpp"

#include <functional>
#include <span>
#include <vector>

namespace nova::fem {

struct DiffusionResidualIndicators2D {
    // Nonnegative elementwise residual estimates.
    std::vector<double> element_errors;
    // Euclidean aggregate: sqrt(sum_K eta_K^2).
    double total_error{};
};

// Residual estimator for -div(k grad(u)) = f on affine triangles.
// Conductivity is positive and constant within each element. Supports P1 and
// P2 H1 fields. It includes the cell residual and interior normal-flux jumps;
// natural-boundary residual terms and coefficient-gradient terms are omitted.
// This is an estimator for refinement, not a guaranteed global error bound.
[[nodiscard]] DiffusionResidualIndicators2D
estimate_diffusion_residual_error_2d(
    const TriangleMesh2D& mesh,
    H1TriangleOrder order,
    std::span<const double> dof_values,
    std::span<const double> element_conductivity,
    const std::function<double(double, double)>& source
);

} // namespace nova::fem
