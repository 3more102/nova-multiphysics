#include "test_framework.hpp"

#include "nova/linalg/bicgstab.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

NOVA_TEST(bicgstab_solves_nonsymmetric_system) {
    const std::vector<nova::linalg::Triplet> entries{
        {0, 0, 4.0}, {0, 1, 1.0},
        {1, 0, -2.0}, {1, 1, 3.0}
    };
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(2, 2, entries);
    const std::vector<double> rhs{6.0, 4.0};

    nova::linalg::SolverOptions options;
    options.relative_tolerance = 1.0e-12;
    options.absolute_tolerance = 1.0e-14;

    const auto result = nova::linalg::solve_bicgstab(matrix, rhs, options);

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE_NEAR(result.solution[0], 1.0, 1.0e-10);
    NOVA_REQUIRE_NEAR(result.solution[1], 2.0, 1.0e-10);
    NOVA_REQUIRE(result.residual_norm < 1.0e-10);

    const auto product = matrix.multiply(result.solution);
    const double true_residual = std::hypot(
        rhs[0] - product[0],
        rhs[1] - product[1]
    );
    NOVA_REQUIRE_NEAR(result.residual_norm, true_residual, 1.0e-14);
}

NOVA_TEST(bicgstab_solves_nonsymmetric_system_with_jacobi_preconditioning) {
    const std::vector<nova::linalg::Triplet> entries{
        {0, 0, 4.0}, {0, 1, 1.0},
        {1, 0, -2.0}, {1, 1, 3.0}
    };
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(2, 2, entries);
    const std::vector<double> rhs{6.0, 4.0};

    nova::linalg::SolverOptions options;
    options.relative_tolerance = 1.0e-12;
    options.absolute_tolerance = 1.0e-14;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;

    const auto result = nova::linalg::solve_bicgstab(matrix, rhs, options);

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE_NEAR(result.solution[0], 1.0, 1.0e-10);
    NOVA_REQUIRE_NEAR(result.solution[1], 2.0, 1.0e-10);
    NOVA_REQUIRE(result.message.find("Jacobi") != std::string::npos);
}

NOVA_TEST(bicgstab_rejects_zero_jacobi_diagonal) {
    const std::vector<nova::linalg::Triplet> entries{
        {0, 1, 1.0},
        {1, 0, 1.0}
    };
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(2, 2, entries);
    const std::vector<double> rhs{1.0, 1.0};

    nova::linalg::SolverOptions options;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;

    bool threw = false;
    try {
        static_cast<void>(nova::linalg::solve_bicgstab(matrix, rhs, options));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}

NOVA_TEST(bicgstab_reports_iteration_limit_without_claiming_convergence) {
    const std::vector<nova::linalg::Triplet> entries{
        {0, 0, 4.0}, {0, 1, 1.0},
        {1, 0, -2.0}, {1, 1, 3.0}
    };
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(2, 2, entries);
    const std::vector<double> rhs{6.0, 4.0};

    nova::linalg::SolverOptions options;
    options.max_iterations = 1;
    options.relative_tolerance = 1.0e-14;
    options.absolute_tolerance = 1.0e-16;

    const auto result = nova::linalg::solve_bicgstab(matrix, rhs, options);

    NOVA_REQUIRE(!result.converged);
    NOVA_REQUIRE(result.iterations == 1);
    NOVA_REQUIRE(result.message.find("maximum iteration") != std::string::npos);
}

NOVA_TEST(bicgstab_rejects_dimension_mismatch) {
    const std::vector<nova::linalg::Triplet> entries{{0, 0, 1.0}};
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(1, 1, entries);
    const std::vector<double> rhs{1.0, 2.0};

    bool threw = false;
    try {
        static_cast<void>(nova::linalg::solve_bicgstab(matrix, rhs));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}
