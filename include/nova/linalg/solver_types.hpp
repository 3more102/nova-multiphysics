#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace nova::linalg {

enum class PreconditionerType {
    None,
    Jacobi,
    IncompleteCholesky0 // CG only
};

struct SolverOptions {
    std::size_t max_iterations{1000};
    double relative_tolerance{1.0e-10};
    double absolute_tolerance{1.0e-14};
    PreconditionerType preconditioner{PreconditionerType::None};
};

struct SolverResult {
    std::vector<double> solution;
    bool converged{false};
    std::size_t iterations{0};
    double initial_residual_norm{0.0};
    double residual_norm{0.0};
    std::string message;
};

} // namespace nova::linalg
