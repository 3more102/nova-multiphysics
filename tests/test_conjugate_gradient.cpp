#include "test_framework.hpp"

#include "nova/linalg/conjugate_gradient.hpp"

#include <array>
#include <cmath>
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

NOVA_TEST(conjugate_gradient_does_not_falsely_converge_on_tiny_relative_rhs) {
    using namespace nova::linalg;
    const auto matrix = CsrMatrix::from_triplets(2, 2, {
        {0,0,4.0}, {0,1,1.0}, {1,0,1.0}, {1,1,3.0}});
    const std::array expected{1e-200, 2e-200};
    const auto rhs = matrix.multiply(expected);
    for (const auto preconditioner : std::array{
             PreconditionerType::None, PreconditionerType::Jacobi}) {
        SolverOptions options;
        options.preconditioner = preconditioner;
        options.absolute_tolerance = 0.0; // Require scale-independent relative accuracy.
        options.relative_tolerance = 1e-11;
        const auto result = solve_conjugate_gradient(matrix, rhs, options);
        NOVA_REQUIRE(result.converged && result.iterations > 0);
        NOVA_REQUIRE(result.initial_residual_norm > 1e-200);
        for (std::size_t i = 0; i < 2; ++i) {
            NOVA_REQUIRE(std::abs(result.solution[i] / expected[i] - 1.0) < 1e-9);
        }
        const auto product = matrix.multiply(result.solution);
        const double true_residual = std::hypot(rhs[0]-product[0], rhs[1]-product[1]);
        NOVA_REQUIRE(result.residual_norm <= 1.00001 * true_residual + 1e-220);
        NOVA_REQUIRE(true_residual <= options.relative_tolerance *
            result.initial_residual_norm);
    }
}

NOVA_TEST(conjugate_gradient_supports_tiny_matrix_with_scaled_preconditioner) {
    using namespace nova::linalg;
    const double scale = 1e-200;
    const auto matrix = CsrMatrix::from_triplets(2, 2, {
        {0,0,4*scale}, {0,1,scale}, {1,0,scale}, {1,1,3*scale}});
    const std::array expected{1.0, 2.0};
    const auto rhs = matrix.multiply(expected);
    SolverOptions options;
    options.preconditioner = PreconditionerType::Jacobi;
    options.absolute_tolerance = 0.0;
    options.relative_tolerance = 1e-11;
    const auto result = solve_conjugate_gradient(matrix, rhs, options);
    NOVA_REQUIRE(result.converged && result.iterations > 0);
    for (std::size_t i = 0; i < 2; ++i) {
        NOVA_REQUIRE_NEAR(result.solution[i], expected[i], 1e-9);
    }
    NOVA_REQUIRE(result.residual_norm <= options.relative_tolerance *
        result.initial_residual_norm);
}

NOVA_TEST(conjugate_gradient_norm_does_not_overflow_on_large_finite_rhs) {
    using namespace nova::linalg;
    const auto matrix = CsrMatrix::from_triplets(2, 2, {{0,0,1.0}, {1,1,1.0}});
    const std::array expected{1e200, -1e200};
    SolverOptions options;
    options.relative_tolerance = 1e-12;
    options.absolute_tolerance = 0.0;
    const auto result = solve_conjugate_gradient(matrix, expected, options);
    NOVA_REQUIRE(result.converged && result.iterations == 1);
    NOVA_REQUIRE(std::isfinite(result.initial_residual_norm));
    NOVA_REQUIRE(std::abs(result.initial_residual_norm / 1e200 - std::sqrt(2.0)) < 1e-13);
    NOVA_REQUIRE(result.solution == std::vector<double>(expected.begin(), expected.end()));
}

NOVA_TEST(conjugate_gradient_reports_accurate_norm_for_subnormal_squared_residual) {
    using namespace nova::linalg;
    const auto identity = CsrMatrix::from_triplets(2, 2, {{0,0,1.0}, {1,1,1.0}});
    const std::array rhs{2e-162, 2e-162};
    SolverOptions options;
    options.absolute_tolerance = 0.0;
    const auto result = solve_conjugate_gradient(identity, rhs, options);
    NOVA_REQUIRE(result.converged && result.iterations == 1);
    NOVA_REQUIRE_NEAR(result.initial_residual_norm / std::hypot(rhs[0], rhs[1]), 1.0, 1e-14);
}

NOVA_TEST(conjugate_gradient_jacobi_accepts_subnormal_positive_diagonal) {
    using namespace nova::linalg;
    // The reciprocal of this finite diagonal overflows, but r/d is exactly 2.
    const double diagonal = 1e-310;
    const auto matrix = CsrMatrix::from_triplets(1, 1, {{0, 0, diagonal}});
    const std::array rhs{2.0 * diagonal};
    SolverOptions options;
    options.preconditioner = PreconditionerType::Jacobi;
    options.absolute_tolerance = 0.0;
    options.relative_tolerance = 1e-10;
    const auto result = solve_conjugate_gradient(matrix, rhs, options);
    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.iterations == 1);
    NOVA_REQUIRE_NEAR(result.solution[0], 2.0, 1e-12);
    const auto product = matrix.multiply(result.solution);
    NOVA_REQUIRE(product[0] == rhs[0]);
    NOVA_REQUIRE(result.residual_norm == 0.0);
}
