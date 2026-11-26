#include "nova/linalg/gmres.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

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
    const double squared = dot(vector, vector);
    if (std::isnormal(squared)) {
        return std::sqrt(squared);
    }
    // Squaring finite components can underflow or overflow even when the
    // Euclidean norm itself is representable. Use a scaled fallback.
    double result = 0.0;
    for (const double value : vector) {
        result = std::hypot(result, value);
    }
    return result;
}

[[nodiscard]] bool all_finite(const std::span<const double> values) {
    return std::all_of(values.begin(), values.end(), [](const double value) {
        return std::isfinite(value);
    });
}

[[nodiscard]] bool near_zero(const double value, const double scale) {
    if (!std::isfinite(value) || !std::isfinite(scale) || scale < 0.0) {
        return true;
    }
    return std::abs(value) <=
        32.0 * std::numeric_limits<double>::epsilon() * scale;
}

} // namespace

SolverResult solve_gmres(
    const CsrMatrix& matrix,
    const std::span<const double> right_hand_side,
    const GmresOptions& options,
    const std::span<const double> initial_guess
) {
    const SolverOptions& solver_options = options.solver;

    if (matrix.rows() == 0 || matrix.rows() != matrix.columns()) {
        throw std::invalid_argument("GMRES requires a non-empty square matrix.");
    }
    if (right_hand_side.size() != matrix.rows()) {
        throw std::invalid_argument("GMRES right-hand-side dimension mismatch.");
    }
    if (!initial_guess.empty() && initial_guess.size() != matrix.columns()) {
        throw std::invalid_argument("GMRES initial-guess dimension mismatch.");
    }
    if (solver_options.max_iterations == 0) {
        throw std::invalid_argument("GMRES iteration limit must be positive.");
    }
    if (options.restart == 0) {
        throw std::invalid_argument("GMRES restart length must be positive.");
    }
    if (!std::isfinite(solver_options.relative_tolerance) ||
        !std::isfinite(solver_options.absolute_tolerance) ||
        solver_options.relative_tolerance < 0.0 ||
        solver_options.absolute_tolerance < 0.0 ||
        (solver_options.relative_tolerance == 0.0 &&
         solver_options.absolute_tolerance == 0.0)) {
        throw std::invalid_argument(
            "GMRES tolerances must be finite, non-negative, and not both zero."
        );
    }
    if (!all_finite(right_hand_side) || !all_finite(initial_guess)) {
        throw std::invalid_argument("GMRES input vectors must contain only finite values.");
    }
    if (solver_options.preconditioner != PreconditionerType::None &&
        solver_options.preconditioner != PreconditionerType::Jacobi) {
        throw std::invalid_argument("Unsupported GMRES preconditioner.");
    }

    std::vector<double> inverse_diagonal;
    if (solver_options.preconditioner == PreconditionerType::Jacobi) {
        const std::vector<double> diagonal = matrix.diagonal();
        inverse_diagonal.resize(diagonal.size(), 0.0);

        for (std::size_t index = 0; index < diagonal.size(); ++index) {
            if (!std::isfinite(diagonal[index]) || diagonal[index] == 0.0) {
                throw std::invalid_argument(
                    "GMRES Jacobi preconditioning requires a finite, nonzero diagonal."
                );
            }

            inverse_diagonal[index] = 1.0 / diagonal[index];
            if (!std::isfinite(inverse_diagonal[index])) {
                throw std::invalid_argument(
                    "GMRES Jacobi inverse diagonal is outside the supported numeric range."
                );
            }
        }
    }

    const auto apply_preconditioner = [&](const std::vector<double>& vector) {
        if (solver_options.preconditioner == PreconditionerType::None) {
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
        throw std::overflow_error("Initial GMRES residual is non-finite.");
    }

    SolverResult result;
    result.initial_residual_norm = norm(residual);
    result.residual_norm = result.initial_residual_norm;

    if (!std::isfinite(result.initial_residual_norm)) {
        throw std::overflow_error("Initial GMRES residual norm overflowed.");
    }

    const double convergence_threshold = std::max(
        solver_options.absolute_tolerance,
        solver_options.relative_tolerance * result.initial_residual_norm
    );
    if (!std::isfinite(convergence_threshold)) {
        throw std::invalid_argument("GMRES convergence threshold overflowed.");
    }

    if (result.residual_norm <= convergence_threshold) {
        result.solution = std::move(solution);
        result.converged = true;
        result.message = "Initial guess satisfies the convergence criterion.";
        return result;
    }

    const std::size_t cycle_length = std::min(
        options.restart,
        solver_options.max_iterations
    );

    while (result.iterations < solver_options.max_iterations) {
        const double beta = norm(residual);
        if (!std::isfinite(beta) || beta == 0.0) {
            result.solution = std::move(solution);
            result.message = "GMRES breakdown: restart residual norm is zero or non-finite.";
            return result;
        }

        std::vector<std::vector<double>> basis(
            cycle_length + 1,
            std::vector<double>(matrix.columns(), 0.0)
        );
        std::vector<std::vector<double>> preconditioned_basis(
            cycle_length,
            std::vector<double>(matrix.columns(), 0.0)
        );
        std::vector<std::vector<double>> hessenberg(
            cycle_length + 1,
            std::vector<double>(cycle_length, 0.0)
        );

        std::vector<double> cosine(cycle_length, 0.0);
        std::vector<double> sine(cycle_length, 0.0);
        std::vector<double> transformed_rhs(cycle_length + 1, 0.0);

        for (std::size_t index = 0; index < residual.size(); ++index) {
            basis[0][index] = residual[index] / beta;
        }
        transformed_rhs[0] = beta;

        std::size_t used_iterations = 0;
        bool happy_breakdown = false;
        bool arnoldi_failure = false;

        for (std::size_t column = 0;
             column < cycle_length &&
             result.iterations < solver_options.max_iterations;
             ++column) {
            preconditioned_basis[column] = apply_preconditioner(basis[column]);
            const std::vector<double> matrix_product =
                matrix.multiply(std::span<const double>(preconditioned_basis[column]));

            if (!all_finite(matrix_product)) {
                result.message = "GMRES breakdown: Arnoldi matrix product became non-finite.";
                arnoldi_failure = true;
                break;
            }

            std::vector<double> work = matrix_product;

            // Two-pass modified Gram-Schmidt reduces loss of orthogonality.
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t row = 0; row <= column; ++row) {
                    const double correction = dot(work, basis[row]);
                    hessenberg[row][column] += correction;

                    for (std::size_t index = 0; index < work.size(); ++index) {
                        work[index] -= correction * basis[row][index];
                    }
                }
            }

            const double next_basis_norm = norm(work);
            const double product_norm = norm(matrix_product);

            if (!std::isfinite(next_basis_norm) || !std::isfinite(product_norm)) {
                result.message = "GMRES breakdown: Arnoldi norm became non-finite.";
                arnoldi_failure = true;
                break;
            }

            hessenberg[column + 1][column] = next_basis_norm;
            happy_breakdown = near_zero(next_basis_norm, product_norm);

            if (!happy_breakdown) {
                for (std::size_t index = 0; index < work.size(); ++index) {
                    basis[column + 1][index] = work[index] / next_basis_norm;
                }
            }

            // Apply all previously computed Givens rotations to this column.
            for (std::size_t rotation = 0; rotation < column; ++rotation) {
                const double upper =
                    cosine[rotation] * hessenberg[rotation][column] +
                    sine[rotation] * hessenberg[rotation + 1][column];
                const double lower =
                    -sine[rotation] * hessenberg[rotation][column] +
                    cosine[rotation] * hessenberg[rotation + 1][column];

                hessenberg[rotation][column] = upper;
                hessenberg[rotation + 1][column] = lower;
            }

            const double diagonal = hessenberg[column][column];
            const double subdiagonal = hessenberg[column + 1][column];
            const double rotation_norm = std::hypot(diagonal, subdiagonal);

            if (!std::isfinite(rotation_norm) || rotation_norm == 0.0) {
                result.message = "GMRES breakdown: singular Givens rotation.";
                arnoldi_failure = true;
                break;
            }

            cosine[column] = diagonal / rotation_norm;
            sine[column] = subdiagonal / rotation_norm;
            hessenberg[column][column] = rotation_norm;
            hessenberg[column + 1][column] = 0.0;

            const double rhs_value = transformed_rhs[column];
            transformed_rhs[column] = cosine[column] * rhs_value;
            transformed_rhs[column + 1] = -sine[column] * rhs_value;

            ++result.iterations;
            used_iterations = column + 1;

            if (std::abs(transformed_rhs[column + 1]) <= convergence_threshold ||
                happy_breakdown) {
                break;
            }
        }

        if (used_iterations == 0) {
            result.solution = std::move(solution);
            if (result.message.empty()) {
                result.message = "GMRES breakdown: no Arnoldi iteration was completed.";
            }
            return result;
        }

        // Back-substitute the upper-triangular least-squares system.
        std::vector<double> coefficients(used_iterations, 0.0);
        for (std::size_t reverse = used_iterations; reverse > 0; --reverse) {
            const std::size_t row = reverse - 1;
            double value = transformed_rhs[row];
            double row_scale = 0.0;

            for (std::size_t column = row; column < used_iterations; ++column) {
                row_scale = std::max(
                    row_scale,
                    std::abs(hessenberg[row][column])
                );
                if (column > row) {
                    value -= hessenberg[row][column] * coefficients[column];
                }
            }

            if (near_zero(hessenberg[row][row], row_scale)) {
                result.solution = std::move(solution);
                result.message = "GMRES breakdown: singular Hessenberg least-squares system.";
                return result;
            }

            coefficients[row] = value / hessenberg[row][row];
            if (!std::isfinite(coefficients[row])) {
                result.solution = std::move(solution);
                result.message = "GMRES breakdown: least-squares coefficient became non-finite.";
                return result;
            }
        }

        for (std::size_t column = 0; column < used_iterations; ++column) {
            for (std::size_t index = 0; index < solution.size(); ++index) {
                solution[index] +=
                    coefficients[column] * preconditioned_basis[column][index];
            }
        }

        if (!all_finite(solution)) {
            result.solution = std::move(solution);
            result.message = "GMRES breakdown: solution became non-finite.";
            return result;
        }

        residual = compute_true_residual(solution);
        if (!all_finite(residual)) {
            result.solution = std::move(solution);
            result.message = "GMRES breakdown: recomputed true residual became non-finite.";
            return result;
        }

        result.residual_norm = norm(residual);
        if (!std::isfinite(result.residual_norm)) {
            result.solution = std::move(solution);
            result.message = "GMRES breakdown: true residual norm became non-finite.";
            return result;
        }

        if (result.residual_norm <= convergence_threshold) {
            result.solution = std::move(solution);
            result.converged = true;
            result.message = solver_options.preconditioner == PreconditionerType::Jacobi
                ? "GMRES converged with Jacobi preconditioning."
                : "GMRES converged.";
            return result;
        }

        if (arnoldi_failure || happy_breakdown) {
            result.solution = std::move(solution);
            if (arnoldi_failure) {
                result.message = "GMRES breakdown before satisfying the true-residual criterion.";
            } else {
                result.message = "GMRES happy breakdown did not satisfy the true-residual criterion.";
            }
            return result;
        }
    }

    result.solution = std::move(solution);
    result.message = "GMRES reached the maximum Arnoldi iteration count before convergence.";
    return result;
}

} // namespace nova::linalg
