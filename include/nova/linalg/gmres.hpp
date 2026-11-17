#pragma once

#include "nova/linalg/csr_matrix.hpp"
#include "nova/linalg/solver_types.hpp"

#include <cstddef>
#include <span>

namespace nova::linalg {

struct GmresOptions {
    SolverOptions solver{};
    std::size_t restart{30};
};

// Solves A*x=b using restarted GMRES for general nonsymmetric square systems.
// Jacobi preconditioning is applied on the right and requires a finite,
// nonzero matrix diagonal. max_iterations counts Arnoldi steps across restarts.
[[nodiscard]] SolverResult solve_gmres(
    const CsrMatrix& matrix,
    std::span<const double> right_hand_side,
    const GmresOptions& options = {},
    std::span<const double> initial_guess = {}
);

} // namespace nova::linalg
