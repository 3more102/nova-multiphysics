#include "nova/linalg/conjugate_gradient.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace nova::linalg {
namespace {

[[nodiscard]] double dot(
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    double result = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        result += left[index] * right[index];
    }
    return result;
}

[[nodiscard]] double norm(const std::vector<double>& vector) {
    return std::sqrt(dot(vector, vector));
}

[[nodiscard]] bool all_finite(const std::span<const double> values) {
    return std::all_of(values.begin(), values.end(), [](const double value) {
        return std::isfinite(value);
    });
}

} // namespace

SolverResult solve_conjugate_gradient(
    const CsrMatrix& matrix,
    const std::span<const double> right_hand_side,
    const SolverOptions& options,
    const std::span<const double> initial_guess
) {
    if (matrix.rows() == 0 || matrix.rows() != matrix.columns()) {
        throw std::invalid_argument("Conjugate gradient requires a non-empty square matrix.");
    }
    if (right_hand_side.size() != matrix.rows()) {
        throw std::invalid_argument("Conjugate-gradient right-hand-side dimension mismatch.");
    }
    if (!initial_guess.empty() && initial_guess.size() != matrix.columns()) {
        throw std::invalid_argument("Conjugate-gradient initial-guess dimension mismatch.");
    }
    if (options.max_iterations == 0) {
        throw std::invalid_argument("Conjugate-gradient iteration limit must be positive.");
    }
    if (!std::isfinite(options.relative_tolerance) ||
        !std::isfinite(options.absolute_tolerance) ||
        options.relative_tolerance < 0.0 ||
        options.absolute_tolerance < 0.0 ||
        (options.relative_tolerance == 0.0 && options.absolute_tolerance == 0.0)) {
        throw std::invalid_argument("Conjugate-gradient tolerances must be finite, non-negative, and not both zero.");
    }
    if (!all_finite(right_hand_side) || !all_finite(initial_guess)) {
        throw std::invalid_argument("Conjugate-gradient input vectors must contain only finite values.");
    }

    std::vector<double> solution(matrix.columns(), 0.0);
    if (!initial_guess.empty()) {
        std::copy(initial_guess.begin(), initial_guess.end(), solution.begin());
    }

    const std::vector<double> matrix_times_solution =
        matrix.multiply(std::span<const double>(solution));

    std::vector<double> residual(right_hand_side.size(), 0.0);
    for (std::size_t index = 0; index < residual.size(); ++index) {
        residual[index] = right_hand_side[index] - matrix_times_solution[index];
    }

    std::vector<double> direction = residual;
    double residual_squared = dot(residual, residual);
    if (!std::isfinite(residual_squared)) {
        throw std::overflow_error("Initial conjugate-gradient residual norm overflowed.");
    }

    SolverResult result;
    result.initial_residual_norm = norm(residual);
    result.residual_norm = result.initial_residual_norm;

    const double convergence_threshold = std::max(
        options.absolute_tolerance,
        options.relative_tolerance * result.initial_residual_norm
    );

    if (result.residual_norm <= convergence_threshold) {
        result.solution = std::move(solution);
        result.converged = true;
        result.message = "Initial guess satisfies the convergence criterion.";
        return result;
    }

    for (std::size_t iteration = 0; iteration < options.max_iterations; ++iteration) {
        const std::vector<double> matrix_times_direction =
            matrix.multiply(std::span<const double>(direction));
        const double curvature = dot(direction, matrix_times_direction);

        if (!std::isfinite(curvature) || curvature <= 0.0) {
            result.solution = std::move(solution);
            result.iterations = iteration;
            result.message = "CG breakdown: non-positive or non-finite search-direction curvature.";
            return result;
        }

        const double alpha = residual_squared / curvature;
        if (!std::isfinite(alpha)) {
            result.solution = std::move(solution);
            result.iterations = iteration;
            result.message = "CG breakdown: non-finite step length.";
            return result;
        }

        for (std::size_t index = 0; index < solution.size(); ++index) {
            solution[index] += alpha * direction[index];
            residual[index] -= alpha * matrix_times_direction[index];
        }

        const double next_residual_squared = dot(residual, residual);
        result.iterations = iteration + 1;
        result.residual_norm = std::sqrt(next_residual_squared);

        if (!std::isfinite(result.residual_norm)) {
            result.solution = std::move(solution);
            result.message = "CG breakdown: residual became non-finite.";
            return result;
        }

        if (result.residual_norm <= convergence_threshold) {
            result.solution = std::move(solution);
            result.converged = true;
            result.message = "Converged.";
            return result;
        }

        if (residual_squared <= 0.0 || !std::isfinite(residual_squared)) {
            result.solution = std::move(solution);
            result.message = "CG breakdown: invalid residual norm.";
            return result;
        }

        const double beta = next_residual_squared / residual_squared;
        if (!std::isfinite(beta)) {
            result.solution = std::move(solution);
            result.message = "CG breakdown: non-finite direction update.";
            return result;
        }

        for (std::size_t index = 0; index < direction.size(); ++index) {
            direction[index] = residual[index] + beta * direction[index];
        }
        residual_squared = next_residual_squared;
    }

    result.solution = std::move(solution);
    result.message = "Maximum iteration count reached before convergence.";
    return result;
}

} // namespace nova::linalg
