#include "nova/nonlinear/newton_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::nonlinear {
namespace {

[[nodiscard]] double vector_norm(const std::vector<double>& values) {
    double squared_norm = 0.0;
    for (const double value : values) {
        squared_norm += value * value;
    }
    return std::sqrt(squared_norm);
}

[[nodiscard]] bool all_finite(const std::span<const double> values) {
    return std::all_of(values.begin(), values.end(), [](const double value) {
        return std::isfinite(value);
    });
}

void validate_residual_shape(
    const std::vector<double>& residual,
    const std::size_t expected_size
) {
    if (residual.size() != expected_size) {
        throw std::invalid_argument("Nonlinear residual dimension does not match the unknown vector.");
    }
}

} // namespace

NewtonResult solve_newton(
    const NonlinearResidualFunction& residual_function,
    const NonlinearJacobianFunction& jacobian_function,
    const std::span<const double> initial_guess,
    const NewtonOptions& options
) {
    if (!residual_function) {
        throw std::invalid_argument("A nonlinear residual function must be provided.");
    }
    if (!jacobian_function) {
        throw std::invalid_argument("A nonlinear Jacobian function must be provided.");
    }
    if (initial_guess.empty()) {
        throw std::invalid_argument("Newton requires at least one unknown.");
    }
    if (!all_finite(initial_guess)) {
        throw std::invalid_argument("Newton initial guess must contain only finite values.");
    }
    if (options.max_iterations == 0) {
        throw std::invalid_argument("Newton iteration limit must be positive.");
    }
    if (!std::isfinite(options.residual_tolerance) ||
        options.residual_tolerance <= 0.0 ||
        !std::isfinite(options.step_tolerance) ||
        options.step_tolerance < 0.0) {
        throw std::invalid_argument("Newton tolerances must be finite; residual tolerance must be positive.");
    }
    if (!std::isfinite(options.armijo_constant) ||
        options.armijo_constant <= 0.0 ||
        options.armijo_constant >= 1.0) {
        throw std::invalid_argument("Armijo constant must be strictly between zero and one.");
    }
    if (!std::isfinite(options.minimum_step_length) ||
        options.minimum_step_length <= 0.0 ||
        options.minimum_step_length > 1.0) {
        throw std::invalid_argument("Minimum Newton step length must be in (0, 1].");
    }
    if (options.linear_solver != NewtonLinearSolver::GMRES &&
        options.linear_solver != NewtonLinearSolver::BiCGSTAB) {
        throw std::invalid_argument("Unsupported Newton linear solver.");
    }

    std::vector<double> solution(initial_guess.begin(), initial_guess.end());
    std::vector<double> residual = residual_function(
        std::span<const double>(solution)
    );
    validate_residual_shape(residual, solution.size());

    if (!all_finite(residual)) {
        throw std::invalid_argument("Initial nonlinear residual contains a non-finite value.");
    }

    NewtonResult result;
    result.solution = solution;
    result.initial_residual_norm = vector_norm(residual);
    result.residual_norm = result.initial_residual_norm;

    if (!std::isfinite(result.initial_residual_norm)) {
        throw std::overflow_error("Initial nonlinear residual norm overflowed.");
    }

    if (result.residual_norm <= options.residual_tolerance) {
        result.converged = true;
        result.message = "Initial guess satisfies the nonlinear residual criterion.";
        return result;
    }

    for (std::size_t iteration = 0; iteration < options.max_iterations; ++iteration) {
        const linalg::CsrMatrix jacobian = jacobian_function(
            std::span<const double>(solution)
        );

        if (jacobian.rows() != solution.size() ||
            jacobian.columns() != solution.size()) {
            throw std::invalid_argument("Nonlinear Jacobian dimensions do not match the unknown vector.");
        }

        std::vector<double> right_hand_side(residual.size(), 0.0);
        for (std::size_t index = 0; index < residual.size(); ++index) {
            right_hand_side[index] = -residual[index];
        }

        if (options.linear_solver == NewtonLinearSolver::GMRES) {
            result.last_linear_solve = linalg::solve_gmres(
                jacobian,
                std::span<const double>(right_hand_side),
                options.gmres
            );
        } else {
            result.last_linear_solve = linalg::solve_bicgstab(
                jacobian,
                std::span<const double>(right_hand_side),
                options.bicgstab
            );
        }

        if (!result.last_linear_solve.converged) {
            result.solution = std::move(solution);
            result.message =
                "Newton linear solve failed: " + result.last_linear_solve.message;
            return result;
        }

        const std::vector<double>& update = result.last_linear_solve.solution;
        if (update.size() != solution.size() || !all_finite(update)) {
            result.solution = std::move(solution);
            result.message = "Newton linear solver returned an invalid update vector.";
            return result;
        }

        const double update_norm = vector_norm(update);
        if (!std::isfinite(update_norm)) {
            result.solution = std::move(solution);
            result.message = "Newton update norm became non-finite.";
            return result;
        }

        bool accepted = false;
        double step_length = 1.0;
        std::vector<double> candidate(solution.size(), 0.0);
        std::vector<double> candidate_residual;
        double candidate_residual_norm = std::numeric_limits<double>::infinity();

        for (std::size_t backtrack = 0;
             backtrack <= options.max_backtracks;
             ++backtrack) {
            for (std::size_t index = 0; index < solution.size(); ++index) {
                candidate[index] = solution[index] + step_length * update[index];
            }

            if (all_finite(candidate)) {
                candidate_residual = residual_function(
                    std::span<const double>(candidate)
                );
                validate_residual_shape(candidate_residual, solution.size());

                if (all_finite(candidate_residual)) {
                    candidate_residual_norm = vector_norm(candidate_residual);

                    if (std::isfinite(candidate_residual_norm) &&
                        candidate_residual_norm <=
                            (1.0 - options.armijo_constant * step_length) *
                                result.residual_norm) {
                        accepted = true;
                        break;
                    }
                }
            }

            step_length *= 0.5;
            if (step_length < options.minimum_step_length) {
                break;
            }
        }

        if (!accepted) {
            result.solution = std::move(solution);
            result.message = "Newton line search failed to find a finite trial with sufficient residual decrease.";
            return result;
        }

        result.last_step_norm = step_length * update_norm;
        solution = std::move(candidate);
        residual = std::move(candidate_residual);
        result.residual_norm = candidate_residual_norm;
        result.iterations = iteration + 1;
        result.solution = solution;

        if (result.residual_norm <= options.residual_tolerance) {
            result.converged = true;
            result.message = "Newton converged.";
            return result;
        }

        const double solution_scale = 1.0 + vector_norm(solution);
        if (result.last_step_norm <= options.step_tolerance * solution_scale) {
            result.message =
                "Newton stagnated: step tolerance reached before residual convergence.";
            return result;
        }
    }

    result.message = "Newton reached the maximum nonlinear iteration count.";
    return result;
}

} // namespace nova::nonlinear
