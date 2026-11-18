#pragma once

#include "nova/fem/mesh2d.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace nova::fem {

// Residual-based element indicators for continuous P1 scalar diffusion with
// positive, elementwise-constant isotropic diffusivity and a volume source.
// The estimator includes cell residuals and interior-edge flux jumps.
// Natural-boundary residuals are intentionally not included.
struct DiffusionResidualEstimate2D {
    std::vector<double> element_indicators;
    double global_estimate{};
};

// Computes the classical residual indicator
//   eta_K^2 = h_K^2 ||f||_L2(K)^2
//           + 1/2 sum_{e interior to K} h_e ||[k grad(u_h) . n_e]||_L2(e)^2.
// For P1 fields and elementwise-constant k, div(k grad(u_h)) is zero inside
// each element. A three-point triangle rule integrates f^2 exactly when f is
// linear or lower order. This is an estimator, not a certified error bound.
[[nodiscard]] DiffusionResidualEstimate2D estimate_p1_diffusion_error_2d(
    const TriangleMesh2D& mesh,
    std::span<const double> nodal_values,
    std::span<const double> element_diffusivity,
    const std::function<double(double, double)>& source
);

// Dörfler bulk marking. Returns element indices in descending indicator
// order (stable by index for ties) until the selected squared indicators
// reach theta times the global squared indicator. Returns empty for a
// zero estimator.
[[nodiscard]] std::vector<std::size_t> mark_elements_dorfler_2d(
    std::span<const double> element_indicators,
    double theta
);

} // namespace nova::fem
