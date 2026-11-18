#pragma once

#include "nova/fem/mesh2d.hpp"

#include <functional>
#include <span>
#include <vector>

namespace nova::fem {

// Residual-based element indicators for continuous P1 scalar diffusion with
// positive, elementwise-constant isotropic diffusivity and a volume source.
// The estimator includes cell residuals and interior-edge flux jumps.
// Natural-boundary residuals are intentionally not included.
struct DiffusionResidualIndicators2D {
    std::vector<double> element_indicators;
    double global_estimate{};
};

// Computes the residual indicator
//   eta_K^2 = h_K^2 ||f||_L2(K)^2
//           + 1/2 sum_{e interior to K} h_e ||[k grad(u_h) . n_e]||_L2(e)^2.
// For P1 fields and elementwise-constant k, div(k grad(u_h)) is zero inside
// each element. The three-point triangle rule integrates f^2 exactly when f
// is linear or lower order. This is a refinement indicator, not a certified
// global error bound.
[[nodiscard]] DiffusionResidualIndicators2D estimate_p1_residual_indicators_2d(
    const TriangleMesh2D& mesh,
    std::span<const double> nodal_values,
    std::span<const double> element_diffusivity,
    const std::function<double(double, double)>& source
);

} // namespace nova::fem
