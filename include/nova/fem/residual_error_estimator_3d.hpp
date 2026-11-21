#pragma once

#include "nova/fem/interior_flux_jumps_3d.hpp"

#include <functional>
#include <span>
#include <vector>

namespace nova::fem {

struct DiffusionResidualIndicators3D {
    // Cell-source residual contribution h_K * ||f||_{L2(K)}.
    std::vector<double> cell_residual_errors;
    // Existing interior normal-flux jump report, shared equally by adjacent cells.
    InteriorFluxJumpReport3D flux_jumps;
    // Combined nonnegative element indicators: sqrt(cell^2 + face_jump^2).
    std::vector<double> element_errors;
    // Euclidean aggregate: sqrt(sum_K eta_K^2).
    double total_error{};
};

// P1 residual indicator for -div(k grad(u_h)) = f on tetrahedra.
// Conductivity is positive and constant within each element. It combines the
// cell-source residual with interior normal-flux jumps. Natural-boundary
// residuals, within-cell coefficient gradients and tensor materials are
// omitted; this is a refinement indicator, not a guaranteed error bound.
[[nodiscard]] DiffusionResidualIndicators3D estimate_diffusion_residual_error_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const double> nodal_values,
    std::span<const double> element_conductivity,
    const std::function<double(double, double, double)>& source
);

} // namespace nova::fem
