#include "test_framework.hpp"

#include "nova/linalg/sparse_symmetric_eigen.hpp"
#include "nova/linalg/symmetric_eigen.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

using nova::linalg::CsrMatrix;
using nova::linalg::Triplet;

[[nodiscard]] CsrMatrix identity(const std::size_t dimension) {
    std::vector<Triplet> entries;
    for (std::size_t i = 0; i < dimension; ++i) entries.push_back({i, i, 1.0});
    return CsrMatrix::from_triplets(dimension, dimension, entries);
}

[[nodiscard]] double mass_inner(
    const std::vector<double>& left,
    const CsrMatrix& mass,
    const std::vector<double>& right) {
    const auto product = mass.multiply(right);
    double result = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) result += left[i] * product[i];
    return result;
}

} // namespace

NOVA_TEST(sparse_generalized_eigen_solves_diagonal_pencil_with_repeated_top_eigenvalue) {
    const auto stiffness = CsrMatrix::from_triplets(4, 4, {
        {0,0,10.0}, {1,1,8.0}, {2,2,5.0}, {3,3,2.0}
    });
    const auto mass = CsrMatrix::from_triplets(4, 4, {
        {0,0,2.0}, {1,1,2.0}, {2,2,1.0}, {3,3,1.0}
    });
    nova::linalg::SparseSymmetricEigenOptions options;
    options.max_modes = 3;
    options.max_subspace = 4;
    const auto result = nova::linalg::solve_sparse_symmetric_generalized_eigen(
        stiffness, mass, options);

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.subspace_dimension == 4);
    NOVA_REQUIRE(result.eigenpairs.size() == 3);
    NOVA_REQUIRE_NEAR(result.eigenpairs[0].eigenvalue, 5.0, 1e-10);
    NOVA_REQUIRE_NEAR(result.eigenpairs[1].eigenvalue, 5.0, 1e-10);
    NOVA_REQUIRE_NEAR(result.eigenpairs[2].eigenvalue, 4.0, 1e-10);
    for (const auto& pair : result.eigenpairs) {
        NOVA_REQUIRE(pair.relative_residual < 1e-11);
        NOVA_REQUIRE_NEAR(mass_inner(pair.eigenvector, mass, pair.eigenvector), 1.0, 1e-10);
    }
    NOVA_REQUIRE(std::abs(mass_inner(
        result.eigenpairs[0].eigenvector, mass, result.eigenpairs[1].eigenvector)) < 1e-10);
}

NOVA_TEST(sparse_generalized_eigen_matches_dense_reference_on_coupled_spd_pencil) {
    const auto stiffness = CsrMatrix::from_triplets(5, 5, {
        {0,0,8.0}, {0,1,-1.0},
        {1,0,-1.0}, {1,1,7.0}, {1,2,-1.0},
        {2,1,-1.0}, {2,2,6.0}, {2,3,-0.5},
        {3,2,-0.5}, {3,3,5.0}, {3,4,-0.25},
        {4,3,-0.25}, {4,4,4.0}
    });
    const auto mass = CsrMatrix::from_triplets(5, 5, {
        {0,0,2.0}, {0,1,0.2},
        {1,0,0.2}, {1,1,2.2}, {1,2,0.1},
        {2,1,0.1}, {2,2,1.8}, {2,3,0.1},
        {3,2,0.1}, {3,3,1.7}, {3,4,0.05},
        {4,3,0.05}, {4,4,1.5}
    });
    const auto reference = nova::linalg::solve_symmetric_generalized_eigen(
        stiffness, mass);
    nova::linalg::SparseSymmetricEigenOptions options;
    options.max_modes = 3;
    options.max_subspace = 5;
    const auto result = nova::linalg::solve_sparse_symmetric_generalized_eigen(
        stiffness, mass, options);

    NOVA_REQUIRE(reference.converged);
    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.eigenpairs.size() == 3);
    for (std::size_t i = 0; i < result.eigenpairs.size(); ++i) {
        const auto& expected = reference.eigenpairs[reference.eigenpairs.size() - 1U - i];
        NOVA_REQUIRE_NEAR(result.eigenpairs[i].eigenvalue, expected.eigenvalue, 1e-9);
        NOVA_REQUIRE(result.eigenpairs[i].relative_residual < 1e-10);
    }
}

NOVA_TEST(sparse_generalized_eigen_reports_unconverged_ritz_pairs_without_faking_success) {
    std::vector<Triplet> stiffness_entries;
    std::vector<Triplet> mass_entries;
    constexpr std::size_t dimension = 40;
    for (std::size_t i = 0; i < dimension; ++i) {
        stiffness_entries.push_back({i, i, static_cast<double>(dimension - i)});
        mass_entries.push_back({i, i, 1.0});
    }
    const auto stiffness = CsrMatrix::from_triplets(dimension, dimension, stiffness_entries);
    const auto mass = CsrMatrix::from_triplets(dimension, dimension, mass_entries);
    nova::linalg::SparseSymmetricEigenOptions options;
    options.max_modes = 2;
    options.max_subspace = 2;
    options.relative_tolerance = 1e-12;
    const auto result = nova::linalg::solve_sparse_symmetric_generalized_eigen(
        stiffness, mass, options);

    NOVA_REQUIRE(!result.converged);
    NOVA_REQUIRE(result.eigenpairs.size() == 2);
    NOVA_REQUIRE(result.subspace_dimension == 2);
    NOVA_REQUIRE(result.eigenpairs[0].relative_residual > options.relative_tolerance);
}

NOVA_TEST(sparse_generalized_eigen_handles_dimensions_above_dense_reference_limit) {
    constexpr std::size_t dimension = 600;
    std::vector<Triplet> stiffness_entries;
    std::vector<Triplet> mass_entries;
    stiffness_entries.reserve(dimension);
    mass_entries.reserve(dimension);
    for (std::size_t i = 0; i < dimension; ++i) {
        stiffness_entries.push_back({i, i, static_cast<double>(dimension - i)});
        mass_entries.push_back({i, i, 1.0});
    }
    const auto stiffness = CsrMatrix::from_triplets(dimension, dimension, stiffness_entries);
    const auto mass = CsrMatrix::from_triplets(dimension, dimension, mass_entries);
    nova::linalg::SparseSymmetricEigenOptions options;
    options.max_modes = 2;
    options.max_subspace = 24;
    options.relative_tolerance = 1e-10;
    const auto result = nova::linalg::solve_sparse_symmetric_generalized_eigen(
        stiffness, mass, options);

    NOVA_REQUIRE(result.eigenpairs.size() == 2);
    NOVA_REQUIRE(result.subspace_dimension <= options.max_subspace);
    NOVA_REQUIRE(result.operator_applications == result.subspace_dimension);
    NOVA_REQUIRE(std::isfinite(result.eigenpairs.front().eigenvalue));
    NOVA_REQUIRE(std::isfinite(result.eigenpairs.front().relative_residual));
}

NOVA_TEST(sparse_generalized_eigen_rejects_invalid_dimensions_symmetry_and_options) {
    const auto identity_matrix = identity(3);
    bool rejected = false;
    try {
        (void)nova::linalg::solve_sparse_symmetric_generalized_eigen(
            CsrMatrix::from_triplets(2, 2, {{0,0,1.0}}), identity_matrix);
    } catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    const auto nonsymmetric = CsrMatrix::from_triplets(2, 2, {
        {0,0,2.0}, {0,1,1.0}, {1,1,2.0}
    });
    rejected = false;
    try {
        (void)nova::linalg::solve_sparse_symmetric_generalized_eigen(
            nonsymmetric, identity(2));
    } catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    const auto indefinite_mass = CsrMatrix::from_triplets(2, 2, {
        {0,0,1.0}, {1,1,-1.0}
    });
    rejected = false;
    try {
        (void)nova::linalg::solve_sparse_symmetric_generalized_eigen(
            identity(2), indefinite_mass);
    } catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    auto options = nova::linalg::SparseSymmetricEigenOptions{};
    options.max_subspace = 0;
    rejected = false;
    try {
        (void)nova::linalg::solve_sparse_symmetric_generalized_eigen(
            identity_matrix, identity_matrix, options);
    } catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);
}


NOVA_TEST(sparse_generalized_eigen_explores_complement_of_invariant_start) {
    // The deterministic first starting vector is proportional to (1, 1).
    // Its exact eigenvalue is 1, but (1, -1) has the largest eigenvalue 3.
    // A zero residual for the first vector alone is not sufficient.
    const auto stiffness = CsrMatrix::from_triplets(2, 2, {
        {0,0,2.0}, {0,1,-1.0}, {1,0,-1.0}, {1,1,2.0}
    });
    const auto mass = identity(2);
    nova::linalg::SparseSymmetricEigenOptions options;
    options.max_modes = 1;
    options.max_subspace = 2;
    const auto result = nova::linalg::solve_sparse_symmetric_generalized_eigen(
        stiffness, mass, options);
    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.subspace_dimension == 2);
    NOVA_REQUIRE(result.operator_applications == 2);
    NOVA_REQUIRE(result.eigenpairs.size() == 1);
    NOVA_REQUIRE_NEAR(result.eigenpairs[0].eigenvalue, 3.0, 1e-12);
    NOVA_REQUIRE(result.eigenpairs[0].relative_residual < 1e-12);
    NOVA_REQUIRE_NEAR(mass_inner(result.eigenpairs[0].eigenvector,
        mass, result.eigenpairs[0].eigenvector), 1.0, 1e-12);
}
