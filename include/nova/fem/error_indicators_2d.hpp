#pragma once

#include "nova/fem/transient_heat_2d.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace nova::fem {

struct H1GradientRecoveryIndicators2D {
    // Nonnegative recovered-gradient error estimate for each triangle.
    std::vector<double> element_errors;
    // Euclidean aggregate: sqrt(sum_K eta_K^2).
    double total_error{};
};

// Zienkiewicz-Zhu-style gradient-recovery indicator. Element gradients are
// area-averaged to mesh vertices; the recovered P1 gradient is compared with
// the P1/P2 finite-element gradient in an elementwise energy norm. This is a
// practical marking indicator, not a guaranteed upper bound on true error.
// An omitted conductivity callback gives unit weighting.
[[nodiscard]] H1GradientRecoveryIndicators2D
estimate_h1_gradient_recovery_error_2d(
    const TriangleMesh2D& mesh,
    H1TriangleOrder order,
    std::span<const double> dof_values,
    std::function<double(double, double)> conductivity = {}
);

// Convenience wrapper using the conductivity field from a transient heat
// problem. The supplied temperature vector must match the requested H1 order.
[[nodiscard]] H1GradientRecoveryIndicators2D
estimate_transient_heat_gradient_recovery_2d(
    const TransientHeat2DProblem& problem,
    H1TriangleOrder order,
    std::span<const double> temperature
);

// Dörfler bulk marking: choose the smallest deterministic set of elements
// whose squared indicators sum to at least fraction * total squared error.
// Returns an empty set when all indicators are zero.
[[nodiscard]] std::vector<std::size_t> mark_elements_dorfler_2d(
    std::span<const double> element_errors,
    double fraction
);

} // namespace nova::fem
