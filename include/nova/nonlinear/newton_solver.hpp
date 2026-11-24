#pragma once

#include "nova/linalg/bicgstab.hpp"
#include "nova/linalg/gmres.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <vector>

namespace nova::nonlinear {

enum class NewtonLinearSolver {
    GMRES,
    BiCGSTAB
};

using NonlinearResidualFunction =
    std::function<std::vector<double>(std::span<const double>)>;

using NonlinearJacobianFunction =
    std::function<linalg::CsrMatrix(std::span<const double>)>;

struct NewtonOptions {
    std::size_t max_iterations{50};
    double residual_tolerance{1.0e-10};
    double step_tolerance{1.0e-12};
    double armijo_constant{1.0e-4};
    double minimum_step_length{1.0e-8};
    std::size_t max_backtracks{24};

    NewtonLinearSolver linear_solver{NewtonLinearSolver::GMRES};
    linalg::GmresOptions gmres{};
    linalg::SolverOptions bicgstab{};
    double relative_residual_tolerance{}; // fraction of initial residual, optional
};

struct NewtonResult {
    std::vector<double> solution;
    bool converged{false};
    std::size_t iterations{0};
    double initial_residual_norm{0.0};
    double residual_norm{0.0};
    double last_step_norm{0.0};
    std::string message;
    linalg::SolverResult last_linear_solve;
};

// Solves F(x)=0 with Newton iterations and an Armijo backtracking line search.
// The Jacobian callback must return the sparse matrix dF/dx at the supplied x.
[[nodiscard]] NewtonResult solve_newton(
    const NonlinearResidualFunction& residual_function,
    const NonlinearJacobianFunction& jacobian_function,
    std::span<const double> initial_guess,
    const NewtonOptions& options = {}
);

} // namespace nova::nonlinear
