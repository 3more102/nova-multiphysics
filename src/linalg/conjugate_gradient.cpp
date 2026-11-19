#include "nova/linalg/conjugate_gradient.hpp"
#include "nova/linalg/incomplete_cholesky.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>

namespace nova::linalg {
namespace {

// A mantissa/exponent inner product prevents tiny Krylov products from
// underflowing when their ratio (alpha or beta) is representable.
struct ScaledDot { double mantissa{}; int exponent{}; };

[[nodiscard]] ScaledDot dot(const std::vector<double>& left, const std::vector<double>& right) {
    double direct = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) direct += left[i] * right[i];
    if (std::isnormal(direct)) {
        int exponent = 0;
        return {std::frexp(direct, &exponent), exponent};
    }
    double left_scale = 0.0, right_scale = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        left_scale = std::max(left_scale, std::abs(left[i]));
        right_scale = std::max(right_scale, std::abs(right[i]));
    }
    if (left_scale == 0.0 || right_scale == 0.0) return {};
    double normalized = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        normalized += (left[i] / left_scale) * (right[i] / right_scale);
    }
    if (normalized == 0.0 || !std::isfinite(normalized)) return {normalized, 0};
    int left_exp = 0, right_exp = 0, sum_exp = 0;
    const double lm = std::frexp(left_scale, &left_exp);
    const double rm = std::frexp(right_scale, &right_exp);
    const double sm = std::frexp(normalized, &sum_exp);
    return {lm * rm * sm, left_exp + right_exp + sum_exp};
}

[[nodiscard]] double dot_ratio(const ScaledDot& numerator, const ScaledDot& denominator) {
    return std::scalbn(numerator.mantissa / denominator.mantissa,
                       numerator.exponent - denominator.exponent);
}

[[nodiscard]] double norm(const std::vector<double>& vector) {
    double squared = 0.0;
    for (double value : vector) squared += value * value;
    if (std::isnormal(squared)) return std::sqrt(squared);
    // Scaled sum of squares (Blue): avoids both squaring underflow and overflow.
    double scale = 0.0, sum = 1.0;
    for (double value : vector) {
        const double magnitude = std::abs(value);
        if (magnitude == 0.0) continue;
        if (magnitude > scale) {
            const double ratio = scale / magnitude;
            sum = 1.0 + sum * ratio * ratio;
            scale = magnitude;
        } else {
            const double ratio = magnitude / scale;
            sum += ratio * ratio;
        }
    }
    return scale * std::sqrt(sum);
}

[[nodiscard]] bool all_finite(const std::span<const double> values) {
    return std::all_of(values.begin(), values.end(), [](const double value) {
        return std::isfinite(value);
    });
}

} // namespace

SolverResult solve_conjugate_gradient_impl(
    const CsrMatrix& matrix,
    const std::span<const double> right_hand_side,
    const SolverOptions& options,
    const std::span<const double> initial_guess,
    const IncompleteCholesky0* reused_factor
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
        throw std::invalid_argument(
            "Conjugate-gradient tolerances must be finite, non-negative, and not both zero."
        );
    }
    if (!all_finite(right_hand_side) || !all_finite(initial_guess)) {
        throw std::invalid_argument("Conjugate-gradient input vectors must contain only finite values.");
    }
    if (options.preconditioner != PreconditionerType::None &&
        options.preconditioner != PreconditionerType::Jacobi &&
        options.preconditioner != PreconditionerType::IncompleteCholesky0) {
        throw std::invalid_argument("Unsupported conjugate-gradient preconditioner.");
    }

    std::optional<IncompleteCholesky0> incomplete_factor;
    std::vector<double> inverse_diagonal;
    if (options.preconditioner == PreconditionerType::Jacobi) {
        const std::vector<double> diagonal = matrix.diagonal();
        inverse_diagonal.resize(diagonal.size(), 0.0);

        for (std::size_t index = 0; index < diagonal.size(); ++index) {
            if (!std::isfinite(diagonal[index]) || diagonal[index] <= 0.0) {
                throw std::invalid_argument(
                    "Jacobi preconditioning requires a finite, strictly positive diagonal."
                );
            }

            inverse_diagonal[index] = 1.0 / diagonal[index];
            if (!std::isfinite(inverse_diagonal[index])) {
                throw std::invalid_argument(
                    "Jacobi inverse diagonal is outside the supported numeric range."
                );
            }
        }
    } else if (options.preconditioner == PreconditionerType::IncompleteCholesky0 &&
               reused_factor == nullptr) {
        incomplete_factor.emplace(matrix);
    } else if (reused_factor != nullptr &&
               (options.preconditioner != PreconditionerType::IncompleteCholesky0 ||
                reused_factor->dimension() != matrix.rows())) {
        throw std::invalid_argument("Cached IC(0) factor does not match the CG system.");
    }

    const auto apply_preconditioner = [&](const std::vector<double>& residual) {
        if (options.preconditioner == PreconditionerType::None) {
            return residual;
        }
        if (reused_factor != nullptr) {
            return reused_factor->apply(residual);
        }
        if (incomplete_factor) {
            return incomplete_factor->apply(residual);
        }

        std::vector<double> result(residual.size(), 0.0);
        for (std::size_t index = 0; index < residual.size(); ++index) {
            result[index] = inverse_diagonal[index] * residual[index];
        }
        return result;
    };

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
    if (!all_finite(residual)) {
        throw std::overflow_error("Initial conjugate-gradient residual is non-finite.");
    }

    SolverResult result;
    result.initial_residual_norm = norm(residual);
    result.residual_norm = result.initial_residual_norm;

    if (!std::isfinite(result.initial_residual_norm)) {
        throw std::overflow_error("Initial conjugate-gradient residual norm overflowed.");
    }

    const double convergence_threshold = std::max(
        options.absolute_tolerance,
        options.relative_tolerance * result.initial_residual_norm
    );
    if (!std::isfinite(convergence_threshold)) {
        throw std::invalid_argument("Conjugate-gradient convergence threshold overflowed.");
    }

    if (result.residual_norm <= convergence_threshold) {
        result.solution = std::move(solution);
        result.converged = true;
        result.message = "Initial guess satisfies the convergence criterion.";
        return result;
    }

    std::vector<double> preconditioned_residual = apply_preconditioner(residual);
    ScaledDot rho = dot(residual, preconditioned_residual);

    if (!std::isfinite(rho.mantissa) || rho.mantissa <= 0.0) {
        result.solution = std::move(solution);
        result.message = "CG breakdown: preconditioned residual product is non-positive or non-finite.";
        return result;
    }

    std::vector<double> direction = preconditioned_residual;

    for (std::size_t iteration = 0; iteration < options.max_iterations; ++iteration) {
        const std::vector<double> matrix_times_direction =
            matrix.multiply(std::span<const double>(direction));
        const ScaledDot curvature = dot(direction, matrix_times_direction);

        if (!std::isfinite(curvature.mantissa) || curvature.mantissa <= 0.0) {
            result.solution = std::move(solution);
            result.iterations = iteration;
            result.message = "CG breakdown: non-positive or non-finite search-direction curvature.";
            return result;
        }

        const double alpha = dot_ratio(rho, curvature);
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

        if (!all_finite(solution) || !all_finite(residual)) {
            result.solution = std::move(solution);
            result.iterations = iteration + 1;
            result.message = "CG breakdown: solution or residual became non-finite.";
            return result;
        }

        result.iterations = iteration + 1;
        result.residual_norm = norm(residual);

        if (!std::isfinite(result.residual_norm)) {
            result.solution = std::move(solution);
            result.message = "CG breakdown: residual norm became non-finite.";
            return result;
        }

        if (result.residual_norm <= convergence_threshold) {
            // Recompute the true residual before declaring convergence. This
            // prevents accumulated recurrence round-off from causing a false pass.
            const std::vector<double> true_matrix_times_solution =
                matrix.multiply(std::span<const double>(solution));
            std::vector<double> true_residual(right_hand_side.size(), 0.0);

            for (std::size_t index = 0; index < true_residual.size(); ++index) {
                true_residual[index] =
                    right_hand_side[index] - true_matrix_times_solution[index];
            }

            if (!all_finite(true_residual)) {
                result.solution = std::move(solution);
                result.message = "CG breakdown: recomputed true residual became non-finite.";
                return result;
            }

            result.residual_norm = norm(true_residual);
            if (!std::isfinite(result.residual_norm)) {
                result.solution = std::move(solution);
                result.message = "CG breakdown: true residual norm became non-finite.";
                return result;
            }

            if (result.residual_norm <= convergence_threshold) {
                result.solution = std::move(solution);
                result.converged = true;
                if (options.preconditioner == PreconditionerType::IncompleteCholesky0) {
                    result.message = "Converged with IC(0) preconditioning.";
                } else if (options.preconditioner == PreconditionerType::Jacobi) {
                    result.message = "Converged with Jacobi preconditioning.";
                } else {
                    result.message = "Converged.";
                }
                return result;
            }

            // The recursive residual drifted below tolerance before the true
            // residual did. Restart from the recomputed residual.
            residual = std::move(true_residual);
            preconditioned_residual = apply_preconditioner(residual);
            rho = dot(residual, preconditioned_residual);

            if (!std::isfinite(rho.mantissa) || rho.mantissa <= 0.0) {
                result.solution = std::move(solution);
                result.message = "CG breakdown after true-residual restart.";
                return result;
            }

            direction = preconditioned_residual;
            continue;
        }

        preconditioned_residual = apply_preconditioner(residual);
        const ScaledDot next_rho = dot(residual, preconditioned_residual);

        if (!std::isfinite(next_rho.mantissa) || next_rho.mantissa <= 0.0) {
            result.solution = std::move(solution);
            result.message = "CG breakdown: preconditioned residual product is non-positive or non-finite.";
            return result;
        }

        const double beta = dot_ratio(next_rho, rho);
        if (!std::isfinite(beta)) {
            result.solution = std::move(solution);
            result.message = "CG breakdown: non-finite direction update.";
            return result;
        }

        for (std::size_t index = 0; index < direction.size(); ++index) {
            direction[index] = preconditioned_residual[index] + beta * direction[index];
        }
        rho = next_rho;
    }

    result.solution = std::move(solution);
    result.message = "Maximum iteration count reached before convergence.";
    return result;
}

SolverResult solve_conjugate_gradient(
    const CsrMatrix& matrix,
    const std::span<const double> right_hand_side,
    const SolverOptions& options,
    const std::span<const double> initial_guess
) {
    return solve_conjugate_gradient_impl(
        matrix, right_hand_side, options, initial_guess, nullptr);
}

SolverResult CachedConjugateGradient::solve(
    const CsrMatrix& matrix,
    const std::span<const double> right_hand_side,
    const SolverOptions& options,
    const std::span<const double> initial_guess
) {
    if (options.preconditioner != PreconditionerType::IncompleteCholesky0) {
        return solve_conjugate_gradient(matrix, right_hand_side, options, initial_guess);
    }

    bool identical = cached_matrix_.has_value() && cached_factor_.has_value() &&
        cached_matrix_->rows() == matrix.rows() &&
        cached_matrix_->columns() == matrix.columns() &&
        cached_matrix_->non_zeros() == matrix.non_zeros();
    if (identical) {
        matrix.for_each_nonzero([&](std::size_t i, std::size_t j, double value) {
            if (identical && cached_matrix_->at(i, j) != value) identical = false;
        });
    }

    if (!identical) {
        // Build everything before replacing the previous valid cache.
        IncompleteCholesky0 next_factor(matrix);
        CsrMatrix next_matrix = matrix;
        cached_matrix_ = std::move(next_matrix);
        cached_factor_ = std::move(next_factor);
        ++factor_builds_;
    }

    return solve_conjugate_gradient_impl(
        matrix, right_hand_side, options, initial_guess, &*cached_factor_);
}

} // namespace nova::linalg
