#pragma once

#include "nova/fem/transient_heat_3d.hpp"

#include <vector>

namespace nova::fem::detail {

struct PreparedTransientHeat3D {
    ScalarH1Operators operators;
    std::vector<double> jacobians;
    std::vector<bool> constrained;
};

// Validate static topology and boundary assignments, then assemble material
// operators exactly once for a complete fixed/adaptive study.
[[nodiscard]] PreparedTransientHeat3D prepare_transient_heat_3d(
    const TransientHeat3DProblem& problem
);

[[nodiscard]] std::vector<double> initial_transient_heat_state_3d(
    const TransientHeat3DProblem& problem
);

// Runs the shared theta integrator with preassembled static operators.
[[nodiscard]] TransientHeat3DResult solve_prepared_transient_heat_3d(
    const TransientHeat3DProblem& problem,
    const TransientHeat3DOptions& options,
    const PreparedTransientHeat3D& prepared,
    std::vector<double> initial_state
);

} // namespace nova::fem::detail
