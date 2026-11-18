#pragma once

#include "nova/fem/transient_heat_2d.hpp"

namespace nova::fem::detail {

// Validate the public fixed-step input contract before assembling expensive operators.
void validate_transient_heat_2d_inputs(
    const TransientHeat2DProblem& problem,
    const TransientHeat2DOptions& options
);

// Internal execution path for callers that have already validated inputs and
// assembled operators for this exact problem mesh, order and material model.
[[nodiscard]] TransientHeat2DResult solve_transient_heat_2d_with_operators(
    const TransientHeat2DProblem& problem,
    const TransientHeat2DOptions& options,
    const ScalarH1Operators2D& operators
);

} // namespace nova::fem::detail
