#include "test_framework.hpp"

#include "nova/linalg/conjugate_gradient.hpp"

#include <stdexcept>
#include <string>
#include <vector>

NOVA_TEST(conjugate_gradient_solves_small_spd_system) {
    const std::vector<nova::linalg::Triplet> entries{
        {0, 0, 4.0}, {0, 1, 1.0},
        {1, 0, 1.0}, {1, 1, 3.0}
    };
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(2, 2, entries);
    const std::vector<double> rhs{1.0, 2.0};

    nova::linalg::SolverOptions options;
    options.relative_tolerance = 1.0e-12;
    options.absolute_tolerance = 1.0e-14;

    const auto result = nova::linalg::solve_conjugate_gradient(matrix, rhs, options);

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.iterations > 0);
    NOVA_REQUIRE_NEAR(result.solution[0], 1.0 / 11.0, 1.0e-11);
    NOVA_REQUIRE_NEAR(result.solution[1], 7.0 / 11.0, 1.0e-11);
    NOVA_REQUIRE(result.residual_norm < 1.0e-10);

    const auto product = matrix.multiply(result.solution);
    const double true_residual =
        std::hypot(rhs[0] - product[0], rhs[1] - product[1]);
    NOVA_REQUIRE_NEAR(result.residual_norm, true_residual, 1.0e-14);
}

NOVA_TEST(conjugate_gradient_jacobi_preconditioner_solves_scaled_diagonal_system) {
    const std::vector<nova::linalg::Triplet> entries{
        {0, 0, 1.0},
        {1, 1, 100.0},
        {2, 2, 10000.0}
    };
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(3, 3, entries);
    const std::vector<double> rhs{1.0, 100.0, 10000.0};

    nova::linalg::SolverOptions options;
    options.relative_tolerance = 1.0e-12;
    options.absolute_tolerance = 1.0e-14;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;

    const auto result = nova::linalg::solve_conjugate_gradient(matrix, rhs, options);

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.iterations == 1);
    NOVA_REQUIRE(result.message.find("Jacobi") != std::string::npos);
    for (const double value : result.solution) {
        NOVA_REQUIRE_NEAR(value, 1.0, 1.0e-12);
    }
    NOVA_REQUIRE(result.residual_norm < 1.0e-12);
}

NOVA_TEST(conjugate_gradient_jacobi_preconditioner_rejects_non_positive_diagonal) {
    const std::vector<nova::linalg::Triplet> entries{{0, 0, 0.0}};
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(1, 1, entries);
    const std::vector<double> rhs{1.0};

    nova::linalg::SolverOptions options;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;

    bool threw = false;
    try {
        static_cast<void>(nova::linalg::solve_conjugate_gradient(matrix, rhs, options));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}

NOVA_TEST(conjugate_gradient_reports_non_positive_curvature) {
    const std::vector<nova::linalg::Triplet> entries{{0, 0, -1.0}};
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(1, 1, entries);
    const std::vector<double> rhs{1.0};

    const auto result = nova::linalg::solve_conjugate_gradient(matrix, rhs);

    NOVA_REQUIRE(!result.converged);
    NOVA_REQUIRE(result.message.find("breakdown") != std::string::npos);
}

NOVA_TEST(conjugate_gradient_rejects_dimension_mismatch) {
    const std::vector<nova::linalg::Triplet> entries{{0, 0, 1.0}};
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(1, 1, entries);
    const std::vector<double> rhs{1.0, 2.0};

    bool threw = false;
    try {
        static_cast<void>(nova::linalg::solve_conjugate_gradient(matrix, rhs));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}
