#pragma once

#include "nova/linalg/csr_matrix.hpp"
#include "nova/linalg/solver_types.hpp"

#include <optional>

namespace nova::fem::detail {

// Consumes the assembled triplets and RHS. Prescribed identity rows receive
// their exact values in the initial guess, including zero-valued gauge pins.
[[nodiscard]] linalg::SolverResult solve_constrained_diffusion(
    std::vector<linalg::Triplet> entries, std::vector<double> rhs,
    const std::vector<std::optional<double>>& prescribed, const linalg::SolverOptions& options
);

} // namespace nova::fem::detail
