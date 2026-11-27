#include "nova/linalg/sparse_symmetric_eigen.hpp"

#include "nova/linalg/conjugate_gradient.hpp"
#include "nova/linalg/symmetric_eigen.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::linalg {
namespace {

[[nodiscard]] double dot(
    const std::vector<double>& left,
    const std::vector<double>& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("Sparse eigen inner-product dimension mismatch.");
    }
    long double sum = 0.0L;
    for (std::size_t i = 0; i < left.size(); ++i) {
        sum += static_cast<long double>(left[i]) *
               static_cast<long double>(right[i]);
    }
    if (!std::isfinite(sum) ||
        std::abs(sum) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Sparse eigen inner product overflow.");
    }
    return static_cast<double>(sum);
}

[[nodiscard]] double vector_norm(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw std::overflow_error("Sparse eigen vector contains a nonfinite value.");
        }
        result = std::hypot(result, value);
    }
    return result;
}

void axpy(
    std::vector<double>& destination,
    const double scale,
    const std::vector<double>& source) {
    if (destination.size() != source.size() || !std::isfinite(scale)) {
        throw std::invalid_argument("Invalid sparse eigen vector update.");
    }
    for (std::size_t i = 0; i < destination.size(); ++i) {
        destination[i] = std::fma(-scale, source[i], destination[i]);
        if (!std::isfinite(destination[i])) {
            throw std::overflow_error("Sparse eigen orthogonalization overflow.");
        }
    }
}

void validate_symmetric(const CsrMatrix& matrix, const char* label) {
    if (matrix.rows() == 0 || matrix.rows() != matrix.columns()) {
        throw std::invalid_argument(std::string(label) + " must be a non-empty square matrix.");
    }
    double scale = 0.0;
    matrix.for_each_nonzero([&](std::size_t, std::size_t, const double value) {
        scale = std::max(scale, std::abs(value));
    });
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * scale;
    bool symmetric = true;
    matrix.for_each_nonzero([&](const std::size_t row, const std::size_t column, const double value) {
        if (std::abs(value - matrix.at(column, row)) > tolerance) {
            symmetric = false;
        }
    });
    if (!symmetric) {
        throw std::invalid_argument(std::string(label) + " must be symmetric.");
    }
}

[[nodiscard]] double frobenius_norm(const CsrMatrix& matrix) {
    double result = 0.0;
    matrix.for_each_nonzero([&](std::size_t, std::size_t, const double value) {
        result = std::hypot(result, value);
    });
    return result;
}


[[nodiscard]] long double quadratic_form_roundoff_scale(
    const CsrMatrix& matrix,
    const std::vector<double>& vector) {
    if (vector.size() != matrix.columns()) {
        throw std::invalid_argument("Sparse eigen quadratic-form dimension mismatch.");
    }
    std::vector<long double> row_absolute_products(vector.size(), 0.0L);
    matrix.for_each_nonzero(
        [&](const std::size_t row, const std::size_t column, const double value) {
            row_absolute_products[row] +=
                std::abs(static_cast<long double>(value)) *
                std::abs(static_cast<long double>(vector[column]));
        });
    long double scale = 0.0L;
    for (std::size_t row = 0; row < vector.size(); ++row) {
        scale += std::abs(static_cast<long double>(vector[row])) *
                 row_absolute_products[row];
    }
    if (!std::isfinite(scale)) {
        throw std::overflow_error("Sparse eigen quadratic-form error scale overflow.");
    }
    return scale;
}

[[nodiscard]] long double quadratic_form_roundoff_tolerance(
    const CsrMatrix& matrix,
    const std::vector<double>& vector) {
    return 256.0L * static_cast<long double>(std::numeric_limits<double>::epsilon()) *
           quadratic_form_roundoff_scale(matrix, vector);
}

[[nodiscard]] std::uint64_t splitmix64(std::uint64_t& state) {
    std::uint64_t value = (state += 0x9e3779b97f4a7c15ULL);
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

[[nodiscard]] std::vector<double> deterministic_seed(
    const std::size_t dimension,
    std::uint64_t& state) {
    std::vector<double> seed(dimension);
    for (double& value : seed) {
        value = (splitmix64(state) & 1ULL) != 0ULL ? 1.0 : -1.0;
    }
    return seed;
}

[[nodiscard]] double mass_norm_squared(
    const std::vector<double>& vector,
    const std::vector<double>& mass_vector) {
    const double value = dot(vector, mass_vector);
    if (!std::isfinite(value)) {
        throw std::overflow_error("Sparse eigen mass norm is nonfinite.");
    }
    return value;
}

} // namespace

SparseSymmetricEigenResult solve_sparse_symmetric_generalized_eigen(
    const CsrMatrix& stiffness,
    const CsrMatrix& mass,
    const SparseSymmetricEigenOptions& options) {
    const std::size_t dimension = stiffness.rows();
    if (dimension == 0U || stiffness.columns() != dimension ||
        mass.rows() != dimension || mass.columns() != dimension) {
        throw std::invalid_argument("Sparse generalized eigen matrices must be non-empty and dimensionally compatible.");
    }
    if (options.max_modes == 0U || options.max_subspace == 0U ||
        options.max_subspace > 512U || options.max_mass_iterations == 0U ||
        !std::isfinite(options.relative_tolerance) ||
        options.relative_tolerance <= 0.0 || options.relative_tolerance >= 1.0 ||
        !std::isfinite(options.mass_solve_relative_tolerance) ||
        options.mass_solve_relative_tolerance <= 0.0 ||
        options.mass_solve_relative_tolerance >= 1.0 ||
        !std::isfinite(options.breakdown_tolerance) ||
        options.breakdown_tolerance <= 0.0 || options.breakdown_tolerance >= 1.0) {
        throw std::invalid_argument("Invalid sparse generalized eigen solver options.");
    }

    validate_symmetric(stiffness, "Sparse generalized eigen stiffness");
    validate_symmetric(mass, "Sparse generalized eigen mass");

    const auto mass_diagonal = mass.diagonal();
    for (const double value : mass_diagonal) {
        if (!(value > 0.0) || !std::isfinite(value)) {
            throw std::invalid_argument("Sparse generalized eigen mass must have a finite positive diagonal.");
        }
    }

    const std::size_t target_modes = std::min(options.max_modes, dimension);
    const std::size_t subspace_limit = std::min(options.max_subspace, dimension);
    if (subspace_limit < target_modes) {
        throw std::invalid_argument("Sparse eigen subspace must be at least as large as the requested mode count.");
    }

    SparseSymmetricEigenResult result;
    const double stiffness_norm = frobenius_norm(stiffness);
    const double mass_norm = frobenius_norm(mass);
    if (!std::isfinite(stiffness_norm) || !std::isfinite(mass_norm) ||
        !(mass_norm > 0.0)) {
        throw std::invalid_argument("Sparse generalized eigen matrices have invalid norms.");
    }

    std::vector<std::vector<double>> basis;
    std::vector<std::vector<double>> mass_basis;
    std::vector<std::vector<double>> operator_basis;
    basis.reserve(subspace_limit);
    mass_basis.reserve(subspace_limit);
    operator_basis.reserve(subspace_limit);

    // Deterministic block starts improve recovery of repeated/clustered eigenvalues.
    std::uint64_t random_state = 0x4e4f56415f4c414eULL;
    for (std::size_t seed_index = 0; seed_index < target_modes; ++seed_index) {
        bool accepted = false;
        for (std::size_t attempt = 0; attempt < 16U * target_modes + 8U; ++attempt) {
            auto candidate = deterministic_seed(dimension, random_state);
            auto mass_candidate = mass.multiply(candidate);
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t previous = 0; previous < basis.size(); ++previous) {
                    const double coefficient = dot(basis[previous], mass_candidate);
                    axpy(candidate, coefficient, basis[previous]);
                    axpy(mass_candidate, coefficient, mass_basis[previous]);
                }
            }
            // Recompute M*q after the two orthogonalization passes. Updating
            // M*q by the same floating-point axpy operations can accumulate a
            // mismatch that is significant when the residual is nearly zero.
            mass_candidate = mass.multiply(candidate);
            const double norm_squared = mass_norm_squared(candidate, mass_candidate);
            const long double roundoff_tolerance =
                quadratic_form_roundoff_tolerance(mass, candidate);
            if (static_cast<long double>(norm_squared) < -roundoff_tolerance) {
                throw std::invalid_argument("Sparse generalized eigen mass is not positive definite.");
            }
            if (!(static_cast<long double>(norm_squared) > roundoff_tolerance)) {
                continue;
            }
            const double inverse_norm = 1.0 / std::sqrt(norm_squared);
            for (double& value : candidate) value *= inverse_norm;
            mass_candidate = mass.multiply(candidate);
            basis.push_back(std::move(candidate));
            mass_basis.push_back(std::move(mass_candidate));
            accepted = true;
            break;
        }
        if (!accepted) {
            throw std::runtime_error("Unable to construct an independent M-orthonormal starting block.");
        }
    }

    SolverOptions mass_options;
    mass_options.max_iterations = options.max_mass_iterations;
    mass_options.relative_tolerance = options.mass_solve_relative_tolerance;
    mass_options.absolute_tolerance = 0.0;
    mass_options.preconditioner = PreconditionerType::Jacobi;

    std::size_t next_coordinate = 0U;

    // Generalized Arnoldi expansion in the M inner product. For a symmetric
    // pencil this is a fully reorthogonalized Lanczos/Rayleigh-Ritz process.
    for (std::size_t current = 0; current < basis.size() &&
         current < subspace_limit; ++current) {
        auto applied = stiffness.multiply(basis[current]);
        const auto solve = solve_conjugate_gradient(mass, applied, mass_options);
        result.mass_solve_iterations += solve.iterations;
        if (!solve.converged) {
            throw std::runtime_error(
                "Sparse generalized eigen mass solve failed: " + solve.message);
        }
        operator_basis.push_back(std::move(applied));

        auto candidate = solve.solution;
        auto mass_candidate = mass.multiply(candidate);
        const double initial_norm_squared =
            mass_norm_squared(candidate, mass_candidate);
        const long double initial_roundoff_tolerance =
            quadratic_form_roundoff_tolerance(mass, candidate);
        if (static_cast<long double>(initial_norm_squared) < -initial_roundoff_tolerance) {
            throw std::runtime_error("Mass solve produced a vector with negative M norm.");
        }
        const double initial_norm =
            std::sqrt(std::max(0.0, initial_norm_squared));

        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t previous = 0; previous < basis.size(); ++previous) {
                const double coefficient = dot(basis[previous], mass_candidate);
                axpy(candidate, coefficient, basis[previous]);
                axpy(mass_candidate, coefficient, mass_basis[previous]);
            }
        }

        // Re-evaluate the mass product from the updated vector. The previous
        // incremental update is algebraically equivalent but can lose
        // positive definiteness through cancellation near Krylov breakdown.
        mass_candidate = mass.multiply(candidate);
        const double residual_norm_squared =
            mass_norm_squared(candidate, mass_candidate);
        const long double residual_roundoff_tolerance =
            quadratic_form_roundoff_tolerance(mass, candidate);
        if (static_cast<long double>(residual_norm_squared) <
            -residual_roundoff_tolerance) {
            throw std::runtime_error("Generalized Krylov orthogonalization lost positive mass norm.");
        }
        const double residual_norm =
            std::sqrt(std::max(0.0, residual_norm_squared));
        const double threshold = options.breakdown_tolerance *
            std::max(initial_norm, std::numeric_limits<double>::min());
        if (residual_norm > threshold && basis.size() < subspace_limit) {
            for (double& value : candidate) value /= residual_norm;
            mass_candidate = mass.multiply(candidate);
            const double corrected_norm_squared =
                mass_norm_squared(candidate, mass_candidate);
            if (!(corrected_norm_squared > 0.0) || !std::isfinite(corrected_norm_squared)) {
                throw std::runtime_error("Generalized Krylov vector normalization failed.");
            }
            const double correction = 1.0 / std::sqrt(corrected_norm_squared);
            for (double& value : candidate) value *= correction;
            mass_candidate = mass.multiply(candidate);
            basis.push_back(std::move(candidate));
            mass_basis.push_back(std::move(mass_candidate));
        }

        // An invariant starting subspace can contain exact Ritz pairs while
        // missing larger eigenvalues. If the Krylov queue is exhausted before
        // the configured limit, explore its M-orthogonal complement.
        if (current + 1U == basis.size() && basis.size() < subspace_limit) {
            bool restarted = false;
            while (next_coordinate < dimension) {
                const std::size_t coordinate = next_coordinate++;
                std::vector<double> restart(dimension, 0.0);
                restart[coordinate] = 1.0;
                auto mass_restart = mass.multiply(restart);
                for (int pass = 0; pass < 2; ++pass) {
                    for (std::size_t previous = 0; previous < basis.size(); ++previous) {
                        const double coefficient = dot(basis[previous], mass_restart);
                        axpy(restart, coefficient, basis[previous]);
                        axpy(mass_restart, coefficient, mass_basis[previous]);
                    }
                }
                mass_restart = mass.multiply(restart);
                const double squared = mass_norm_squared(restart, mass_restart);
                const long double tolerance = quadratic_form_roundoff_tolerance(mass, restart);
                if (static_cast<long double>(squared) < -tolerance) {
                    throw std::runtime_error("Restart vector has negative generalized mass norm.");
                }
                if (!(static_cast<long double>(squared) > tolerance) ||
                    std::sqrt(squared) <= options.breakdown_tolerance *
                        std::sqrt(mass_diagonal[coordinate])) {
                    continue;
                }
                const double norm = std::sqrt(squared);
                for (double& value : restart) value /= norm;
                mass_restart = mass.multiply(restart);
                basis.push_back(std::move(restart));
                mass_basis.push_back(std::move(mass_restart));
                restarted = true;
                break;
            }
            if (!restarted) {
                throw std::runtime_error("Unable to expand the exhausted generalized Krylov subspace.");
            }
        }
    }

    const std::size_t projected_dimension = basis.size();
    if (operator_basis.size() != projected_dimension) {
        throw std::runtime_error("Generalized Krylov basis was not fully processed.");
    }
    result.subspace_dimension = projected_dimension;
    result.operator_applications = operator_basis.size();

    std::vector<Triplet> projected_entries;
    if (projected_dimension >
        std::numeric_limits<std::size_t>::max() / projected_dimension) {
        throw std::length_error("Projected generalized eigen matrix is too large.");
    }
    projected_entries.reserve(projected_dimension * projected_dimension);
    for (std::size_t row = 0; row < projected_dimension; ++row) {
        for (std::size_t column = row; column < projected_dimension; ++column) {
            const double forward = dot(basis[row], operator_basis[column]);
            const double reverse = dot(basis[column], operator_basis[row]);
            const double value = 0.5 * forward + 0.5 * reverse;
            if (!std::isfinite(value)) {
                throw std::overflow_error("Projected generalized eigen matrix is nonfinite.");
            }
            if (value != 0.0) {
                projected_entries.push_back({row, column, value});
                if (row != column) projected_entries.push_back({column, row, value});
            }
        }
    }

    const auto projected_stiffness = CsrMatrix::from_triplets(
        projected_dimension, projected_dimension, projected_entries);
    std::vector<Triplet> identity_entries;
    identity_entries.reserve(projected_dimension);
    for (std::size_t index = 0; index < projected_dimension; ++index) {
        identity_entries.push_back({index, index, 1.0});
    }
    const auto projected_mass = CsrMatrix::from_triplets(
        projected_dimension, projected_dimension, identity_entries);
    SymmetricEigenOptions dense_options;
    dense_options.max_dimension = 512U;
    dense_options.max_sweeps = 100U;
    dense_options.relative_tolerance = 1.0e-13;
    const auto projected = solve_symmetric_generalized_eigen(
        projected_stiffness, projected_mass, dense_options);

    const std::size_t returned_modes = std::min(target_modes, projected.eigenpairs.size());
    result.eigenpairs.reserve(returned_modes);
    bool residuals_converged = projected.converged;
    for (std::size_t mode = 0; mode < returned_modes; ++mode) {
        const auto& ritz = projected.eigenpairs[projected.eigenpairs.size() - 1U - mode];
        SparseSymmetricEigenpair pair;
        pair.eigenvalue = ritz.eigenvalue;
        if (!std::isfinite(pair.eigenvalue)) {
            throw std::runtime_error("Projected generalized eigenvalue is nonfinite.");
        }
        pair.eigenvector.assign(dimension, 0.0);
        for (std::size_t column = 0; column < projected_dimension; ++column) {
            axpy(pair.eigenvector, -ritz.eigenvector[column], basis[column]);
        }
        const auto mass_vector = mass.multiply(pair.eigenvector);
        const double norm_squared = mass_norm_squared(pair.eigenvector, mass_vector);
        if (!(norm_squared > 0.0) || !std::isfinite(norm_squared)) {
            throw std::runtime_error("Ritz vector has invalid generalized mass norm.");
        }
        const double inverse_norm = 1.0 / std::sqrt(norm_squared);
        for (double& value : pair.eigenvector) value *= inverse_norm;

        const auto applied_stiffness = stiffness.multiply(pair.eigenvector);
        const auto normalized_mass_vector = mass.multiply(pair.eigenvector);
        std::vector<double> residual(dimension);
        for (std::size_t index = 0; index < dimension; ++index) {
            residual[index] = applied_stiffness[index] -
                pair.eigenvalue * normalized_mass_vector[index];
        }
        const double denominator =
            (stiffness_norm + std::abs(pair.eigenvalue) * mass_norm) *
            vector_norm(pair.eigenvector);
        const double numerator = vector_norm(residual);
        if (!std::isfinite(denominator) || !std::isfinite(numerator)) {
            throw std::overflow_error("Generalized Ritz residual scale overflow.");
        }
        pair.relative_residual = denominator > 0.0 ? numerator / denominator : numerator;
        if (!std::isfinite(pair.relative_residual)) {
            throw std::overflow_error("Generalized Ritz residual is nonfinite.");
        }
        if (pair.relative_residual > options.relative_tolerance) {
            residuals_converged = false;
        }

        const auto largest = std::max_element(
            pair.eigenvector.begin(), pair.eigenvector.end(),
            [](const double left, const double right) {
                return std::abs(left) < std::abs(right);
            });
        if (largest != pair.eigenvector.end() && *largest < 0.0) {
            for (double& value : pair.eigenvector) value = -value;
        }
        result.eigenpairs.push_back(std::move(pair));
    }

    result.converged = residuals_converged &&
        result.eigenpairs.size() == target_modes;
    return result;
}

} // namespace nova::linalg
