#pragma once

#include "nova/linalg/csr_matrix.hpp"

#include <cstddef>
#include <span>
#include <string>
#include <vector>

namespace nova::linalg {

struct SolverOptions {
    std::size_t max_iterations{1000};
    double relative_tolerance{1.0e-10};
    double absolute_tolerance{1.0e-14};
};

struct SolverResult {
    std::vector<double> solution;
    bool converged{false};
    std::size_t iterations{0};
    double initial_residual_norm{0.0};
    double residual_norm{0.0};
    std::string message;
};

// Solves A*x=b using unpreconditioned conjugate gradients.
// The caller is responsible for supplying a symmetric positive-definite matrix.
[[nodiscard]] SolverResult solve_conjugate_gradient(
    const CsrMatrix& matrix,
    std::span<const double> right_hand_side,
    const SolverOptions& options = {},
    std::span<const double> initial_guess = {}
);

} // namespace nova::linalg
