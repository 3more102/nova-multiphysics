#include "nova/linalg/bicgstab.hpp"
#include "krylov_dot.hpp"
#include "krylov_norm.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace nova::linalg {
namespace {

[[nodiscard]] bool all_finite(const std::span<const double> values) {
    return std::all_of(values.begin(), values.end(), [](const double value) {
        return std::isfinite(value);
    });
}

// Test orthogonality independently of the units of either vector.
[[nodiscard]] bool near_orthogonal(
    const detail::ScaledDot& product, double left_norm, double right_norm
) {
    if (!(left_norm > 0.0) || !(right_norm > 0.0) ||
        !std::isfinite(left_norm) || !std::isfinite(right_norm)) {
        return true;
    }
    const double cosine = detail::dot_ratio(product, detail::norm_product(left_norm, right_norm));
    return !std::isfinite(cosine) ||
           std::abs(cosine) <= 16.0 * std::numeric_limits<double>::epsilon();
}

} // namespace

SolverResult solve_bicgstab(
    const CsrMatrix& matrix,
    const std::span<const double> right_hand_side,
    const SolverOptions& options,
    const std::span<const double> initial_guess
) {
    if (matrix.rows() == 0 || matrix.rows() != matrix.columns()) {
        throw std::invalid_argument("BiCGSTAB requires a non-empty square matrix.");
    }
    if (right_hand_side.size() != matrix.rows()) {
        throw std::invalid_argument("BiCGSTAB right-hand-side dimension mismatch.");
    }
    if (!initial_guess.empty() && initial_guess.size() != matrix.columns()) {
        throw std::invalid_argument("BiCGSTAB initial-guess dimension mismatch.");
    }
    if (options.max_iterations == 0) {
        throw std::invalid_argument("BiCGSTAB iteration limit must be positive.");
    }
    if (!std::isfinite(options.relative_tolerance) ||
        !std::isfinite(options.absolute_tolerance) ||
        options.relative_tolerance < 0.0 ||
        options.absolute_tolerance < 0.0 ||
        (options.relative_tolerance == 0.0 && options.absolute_tolerance == 0.0)) {
        throw std::invalid_argument(
            "BiCGSTAB tolerances must be finite, non-negative, and not both zero."
        );
    }
    if (!all_finite(right_hand_side) || !all_finite(initial_guess)) {
        throw std::invalid_argument("BiCGSTAB input vectors must contain only finite values.");
    }
    if (options.preconditioner != PreconditionerType::None &&
        options.preconditioner != PreconditionerType::Jacobi) {
        throw std::invalid_argument("Unsupported BiCGSTAB preconditioner.");
    }

    std::vector<double> inverse_diagonal;
    if (options.preconditioner == PreconditionerType::Jacobi) {
        const std::vector<double> diagonal = matrix.diagonal();
        inverse_diagonal.resize(diagonal.size(), 0.0);

        for (std::size_t index = 0; index < diagonal.size(); ++index) {
            if (!std::isfinite(diagonal[index]) || diagonal[index] == 0.0) {
                throw std::invalid_argument(
                    "BiCGSTAB Jacobi preconditioning requires a finite, nonzero diagonal."
                );
            }

            inverse_diagonal[index] = 1.0 / diagonal[index];
            if (!std::isfinite(inverse_diagonal[index])) {
                throw std::invalid_argument(
                    "BiCGSTAB Jacobi inverse diagonal is outside the supported numeric range."
                );
            }
        }
    }

    const auto apply_preconditioner = [&](const std::vector<double>& vector) {
        if (options.preconditioner == PreconditionerType::None) {
            return vector;
        }

        std::vector<double> result(vector.size(), 0.0);
        for (std::size_t index = 0; index < vector.size(); ++index) {
            result[index] = inverse_diagonal[index] * vector[index];
        }
        return result;
    };

    const auto compute_true_residual = [&](
        const std::vector<double>& solution
    ) {
        const std::vector<double> product =
            matrix.multiply(std::span<const double>(solution));
        std::vector<double> residual(right_hand_side.size(), 0.0);

        for (std::size_t index = 0; index < residual.size(); ++index) {
            residual[index] = right_hand_side[index] - product[index];
        }
        return residual;
    };

    std::vector<double> solution(matrix.columns(), 0.0);
    if (!initial_guess.empty()) {
        std::copy(initial_guess.begin(), initial_guess.end(), solution.begin());
    }

    std::vector<double> residual = compute_true_residual(solution);
    if (!all_finite(residual)) {
        throw std::overflow_error("Initial BiCGSTAB residual is non-finite.");
    }

    SolverResult result;
    result.initial_residual_norm = detail::euclidean_norm(residual);
    result.residual_norm = result.initial_residual_norm;

    if (!std::isfinite(result.initial_residual_norm)) {
        throw std::overflow_error("Initial BiCGSTAB residual norm overflowed.");
    }

    const double convergence_threshold = std::max(
        options.absolute_tolerance,
        options.relative_tolerance * result.initial_residual_norm
    );
    if (!std::isfinite(convergence_threshold)) {
        throw std::invalid_argument("BiCGSTAB convergence threshold overflowed.");
    }

    if (result.residual_norm <= convergence_threshold) {
        result.solution = std::move(solution);
        result.converged = true;
        result.message = "Initial guess satisfies the convergence criterion.";
        return result;
    }

    std::vector<double> shadow_residual;
    std::vector<double> direction(residual.size(), 0.0);
    std::vector<double> previous_product(residual.size(), 0.0);

    detail::ScaledDot previous_rho{};
    double alpha = 1.0;
    double omega = 1.0;
    double residual_unit = 1.0;
    double shadow_norm = 0.0;
    double residual_norm = 0.0;
    bool restart = true;

    for (std::size_t iteration = 0; iteration < options.max_iterations; ++iteration) {
        if (restart) {
            // Solve for corrections in units of the physical residual. A common
            // rescaling of A and b would otherwise square that scale in A*p.
            residual_unit = std::scalbn(1.0, std::ilogb(result.residual_norm));
            for (double& value : residual) value /= residual_unit;
            shadow_residual = residual;
            residual_norm = shadow_norm = detail::euclidean_norm(residual);
        }
        const auto rho = detail::dot(shadow_residual, residual);

        if (near_orthogonal(rho, shadow_norm, residual_norm)) {
            result.solution = std::move(solution);
            result.iterations = iteration;
            result.message = "BiCGSTAB breakdown: rho became zero or numerically singular.";
            return result;
        }

        if (restart) {
            direction = residual;
            restart = false;
        } else {
            const double beta = detail::dot_ratio(rho, previous_rho) * (alpha / omega);
            if (!std::isfinite(beta)) {
                result.solution = std::move(solution);
                result.iterations = iteration;
                result.message = "BiCGSTAB breakdown: beta became non-finite.";
                return result;
            }

            for (std::size_t index = 0; index < direction.size(); ++index) {
                direction[index] = residual[index] +
                    beta * (direction[index] - omega * previous_product[index]);
            }
        }

        const std::vector<double> preconditioned_direction =
            apply_preconditioner(direction);
        const std::vector<double> product =
            matrix.multiply(std::span<const double>(preconditioned_direction));

        const auto denominator = detail::dot(shadow_residual, product);

        if (near_orthogonal(denominator, shadow_norm, detail::euclidean_norm(product))) {
            result.solution = std::move(solution);
            result.iterations = iteration;
            result.message = "BiCGSTAB breakdown: alpha denominator became zero or numerically singular.";
            return result;
        }

        alpha = detail::dot_ratio(rho, denominator);
        if (!std::isfinite(alpha)) {
            result.solution = std::move(solution);
            result.iterations = iteration;
            result.message = "BiCGSTAB breakdown: alpha became non-finite.";
            return result;
        }

        std::vector<double> intermediate_residual(residual.size(), 0.0);
        for (std::size_t index = 0; index < residual.size(); ++index) {
            intermediate_residual[index] = residual[index] - alpha * product[index];
        }

        if (!all_finite(intermediate_residual)) {
            result.solution = std::move(solution);
            result.iterations = iteration + 1;
            result.message = "BiCGSTAB breakdown: intermediate residual became non-finite.";
            return result;
        }

        const double intermediate_norm = detail::euclidean_norm(intermediate_residual);
        if (!std::isfinite(intermediate_norm)) {
            result.solution = std::move(solution);
            result.iterations = iteration + 1;
            result.message = "BiCGSTAB breakdown: intermediate residual norm became non-finite.";
            return result;
        }

        if (residual_unit * intermediate_norm <= convergence_threshold) {
            for (std::size_t index = 0; index < solution.size(); ++index) {
                solution[index] += residual_unit * (alpha * preconditioned_direction[index]);
            }

            std::vector<double> true_residual = compute_true_residual(solution);
            if (!all_finite(true_residual)) {
                result.solution = std::move(solution);
                result.iterations = iteration + 1;
                result.message = "BiCGSTAB breakdown: recomputed true residual became non-finite.";
                return result;
            }

            result.residual_norm = detail::euclidean_norm(true_residual);
            result.iterations = iteration + 1;

            if (std::isfinite(result.residual_norm) &&
                result.residual_norm <= convergence_threshold) {
                result.solution = std::move(solution);
                result.converged = true;
                result.message = options.preconditioner == PreconditionerType::Jacobi
                    ? "Converged with Jacobi-preconditioned BiCGSTAB."
                    : "BiCGSTAB converged.";
                return result;
            }

            if (!std::isfinite(result.residual_norm)) {
                result.solution = std::move(solution);
                result.message = "BiCGSTAB breakdown: true residual norm became non-finite.";
                return result;
            }

            residual = std::move(true_residual);
            restart = true;
            continue;
        }

        const std::vector<double> preconditioned_intermediate =
            apply_preconditioner(intermediate_residual);
        const std::vector<double> second_product =
            matrix.multiply(std::span<const double>(preconditioned_intermediate));

        const double second_norm = detail::euclidean_norm(second_product);
        const auto stabilization_product = detail::dot(second_product, intermediate_residual);
        if (near_orthogonal(stabilization_product, second_norm, intermediate_norm)) {
            result.solution = std::move(solution);
            result.iterations = iteration + 1;
            result.message = "BiCGSTAB breakdown: stabilization product is zero or numerically orthogonal.";
            return result;
        }

        omega = detail::dot_ratio(stabilization_product, detail::norm_product(second_norm, second_norm));
        if (!std::isfinite(omega) || omega == 0.0) {
            result.solution = std::move(solution);
            result.iterations = iteration + 1;
            result.message = "BiCGSTAB breakdown: omega became zero or non-finite.";
            return result;
        }

        for (std::size_t index = 0; index < solution.size(); ++index) {
            solution[index] += residual_unit * (
                alpha * preconditioned_direction[index] +
                omega * preconditioned_intermediate[index]
            );

            residual[index] =
                intermediate_residual[index] - omega * second_product[index];
        }

        if (!all_finite(solution) || !all_finite(residual)) {
            result.solution = std::move(solution);
            result.iterations = iteration + 1;
            result.message = "BiCGSTAB breakdown: solution or residual became non-finite.";
            return result;
        }

        result.iterations = iteration + 1;
        residual_norm = detail::euclidean_norm(residual);
        result.residual_norm = residual_unit * residual_norm;

        if (!std::isfinite(result.residual_norm)) {
            result.solution = std::move(solution);
            result.message = "BiCGSTAB breakdown: residual norm became non-finite.";
            return result;
        }

        if (result.residual_norm <= convergence_threshold) {
            std::vector<double> true_residual = compute_true_residual(solution);
            if (!all_finite(true_residual)) {
                result.solution = std::move(solution);
                result.message = "BiCGSTAB breakdown: recomputed true residual became non-finite.";
                return result;
            }

            result.residual_norm = detail::euclidean_norm(true_residual);
            if (!std::isfinite(result.residual_norm)) {
                result.solution = std::move(solution);
                result.message = "BiCGSTAB breakdown: true residual norm became non-finite.";
                return result;
            }

            if (result.residual_norm <= convergence_threshold) {
                result.solution = std::move(solution);
                result.converged = true;
                result.message = options.preconditioner == PreconditionerType::Jacobi
                    ? "Converged with Jacobi-preconditioned BiCGSTAB."
                    : "BiCGSTAB converged.";
                return result;
            }

            residual = std::move(true_residual);
            restart = true;
            continue;
        }

        previous_rho = rho;
        previous_product = product;
    }

    result.solution = std::move(solution);
    result.message = "BiCGSTAB reached the maximum iteration count before convergence.";
    return result;
}

} // namespace nova::linalg
