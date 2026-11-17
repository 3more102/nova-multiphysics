#pragma once

#include "nova/linalg/csr_matrix.hpp"
#include "nova/linalg/solver_types.hpp"

#include <span>

namespace nova::linalg {

// Solves A*x=b using BiCGSTAB for general nonsymmetric square systems.
// Jacobi preconditioning requires a finite, nonzero matrix diagonal.
[[nodiscard]] SolverResult solve_bicgstab(
    const CsrMatrix& matrix,
    std::span<const double> right_hand_side,
    const SolverOptions& options = {},
    std::span<const double> initial_guess = {}
);

} // namespace nova::linalg
