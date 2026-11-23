#pragma once

#include "nova/physics/electrothermal_3d.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::physics {

// Experimental stationary two-way electrothermal coupling by relaxed Picard
// iteration. Only cellwise scalar or SPD-tensor reference electrical
// conductivities are supported; sigma_cell(T) = sigma_ref_cell * scale(Tbar).
// Tbar is the arithmetic average of the four P1 nodal temperatures in a cell.
// The scalar positive scale is dimensionless and takes temperature in kelvin.
struct NonlinearElectrothermal3DProblem {
    Electrothermal3DProblem base;
    // Required positive electrical conductivity multiplier sigma(T)/sigma_ref.
    std::function<double(double)> conductivity_scale;
    // Optional positive thermal conductivity multiplier k(T)/k_ref.
    // When empty, the supplied thermal material remains unchanged.
    // When set, thermal material must be a complete per-cell scalar or tensor.
    std::function<double(double)> thermal_conductivity_scale;
    double initial_temperature_kelvin{300.0};
};

struct NonlinearElectrothermal3DOptions {
    std::size_t max_iterations{40};
    double relaxation{0.7};  // (0, 1], applied to the temperature iterate
    double absolute_tolerance_kelvin{1.0e-8};
    double relative_tolerance{1.0e-9};
};

struct NonlinearElectrothermal3DResult {
    // Fields evaluated with the last temperature iterate, not with a
    // post-hoc conductivity approximation.
    Electrothermal3DResult fields;
    // Temperature used to evaluate electrical and (if enabled) thermal
    // conductivity on the last solve.
    std::vector<double> last_temperature_iterate;
    std::size_t iterations{};
    // Difference between raw solved temperature and last iterate, before
    // relaxation (avoids false convergence at very small relaxation).
    double max_fixed_point_residual_kelvin{};
    double max_normalized_residual{};
    bool converged{};
};

// The input is never modified. Nonconvergence returns converged=false with
// final diagnostics and fields; invalid inputs and failed linear solves throw.
// The convergence criterion is the unrelaxed nodal fixed-point mismatch:
// |T_solved-T_iter| <= abs_tol + rel_tol * max(|T_solved|, |T_iter|).
// This is not a certified nonlinear residual/error estimator.
[[nodiscard]] NonlinearElectrothermal3DResult solve_nonlinear_electrothermal_3d(
    const NonlinearElectrothermal3DProblem& problem,
    const NonlinearElectrothermal3DOptions& options = {},
    const linalg::SolverOptions& electrical_options = {},
    const linalg::SolverOptions& thermal_options = {}
);

} // namespace nova::physics
