#include "test_framework.hpp"

#include "nova/linalg/gmres.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

nova::linalg::CsrMatrix make_nonsymmetric_matrix() {
    const std::vector<nova::linalg::Triplet> entries{
        {0, 0, 4.0}, {0, 1, 1.0},
        {1, 0, -2.0}, {1, 1, 3.0}, {1, 2, 1.0},
        {2, 1, 1.0}, {2, 2, 2.0}
    };
    return nova::linalg::CsrMatrix::from_triplets(3, 3, entries);
}

} // namespace

NOVA_TEST(gmres_solves_nonsymmetric_system) {
    const auto matrix = make_nonsymmetric_matrix();
    const std::vector<double> rhs{6.0, 3.0, 0.0};

    nova::linalg::GmresOptions options;
    options.restart = 3;
    options.solver.relative_tolerance = 1.0e-12;
    options.solver.absolute_tolerance = 1.0e-14;

    const auto result = nova::linalg::solve_gmres(matrix, rhs, options);

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE_NEAR(result.solution[0], 1.0, 1.0e-10);
    NOVA_REQUIRE_NEAR(result.solution[1], 2.0, 1.0e-10);
    NOVA_REQUIRE_NEAR(result.solution[2], -1.0, 1.0e-10);
    NOVA_REQUIRE(result.residual_norm < 1.0e-10);

    const auto product = matrix.multiply(result.solution);
    const double true_residual = std::sqrt(
        (rhs[0] - product[0]) * (rhs[0] - product[0]) +
        (rhs[1] - product[1]) * (rhs[1] - product[1]) +
        (rhs[2] - product[2]) * (rhs[2] - product[2])
    );
    NOVA_REQUIRE_NEAR(result.residual_norm, true_residual, 1.0e-14);
}

NOVA_TEST(gmres_restarts_krylov_subspace_and_converges) {
    const auto matrix = make_nonsymmetric_matrix();
    const std::vector<double> rhs{6.0, 3.0, 0.0};

    nova::linalg::GmresOptions options;
    options.restart = 2;
    options.solver.max_iterations = 50;
    options.solver.relative_tolerance = 1.0e-10;
    options.solver.absolute_tolerance = 1.0e-12;

    const auto result = nova::linalg::solve_gmres(matrix, rhs, options);

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.iterations > options.restart);
    NOVA_REQUIRE_NEAR(result.solution[0], 1.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.solution[1], 2.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.solution[2], -1.0, 1.0e-8);
}

NOVA_TEST(gmres_supports_jacobi_preconditioning) {
    const auto matrix = make_nonsymmetric_matrix();
    const std::vector<double> rhs{6.0, 3.0, 0.0};

    nova::linalg::GmresOptions options;
    options.restart = 3;
    options.solver.relative_tolerance = 1.0e-12;
    options.solver.absolute_tolerance = 1.0e-14;
    options.solver.preconditioner = nova::linalg::PreconditionerType::Jacobi;

    const auto result = nova::linalg::solve_gmres(matrix, rhs, options);

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE_NEAR(result.solution[0], 1.0, 1.0e-10);
    NOVA_REQUIRE_NEAR(result.solution[1], 2.0, 1.0e-10);
    NOVA_REQUIRE_NEAR(result.solution[2], -1.0, 1.0e-10);
    NOVA_REQUIRE(result.message.find("Jacobi") != std::string::npos);
}

NOVA_TEST(gmres_rejects_zero_jacobi_diagonal) {
    const std::vector<nova::linalg::Triplet> entries{
        {0, 1, 1.0},
        {1, 0, 1.0}
    };
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(2, 2, entries);
    const std::vector<double> rhs{1.0, 1.0};

    nova::linalg::GmresOptions options;
    options.solver.preconditioner = nova::linalg::PreconditionerType::Jacobi;

    bool threw = false;
    try {
        static_cast<void>(nova::linalg::solve_gmres(matrix, rhs, options));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}

NOVA_TEST(gmres_rejects_zero_restart_length) {
    const auto matrix = make_nonsymmetric_matrix();
    const std::vector<double> rhs{6.0, 3.0, 0.0};

    nova::linalg::GmresOptions options;
    options.restart = 0;

    bool threw = false;
    try {
        static_cast<void>(nova::linalg::solve_gmres(matrix, rhs, options));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}

NOVA_TEST(gmres_reports_iteration_limit_without_claiming_convergence) {
    const auto matrix = make_nonsymmetric_matrix();
    const std::vector<double> rhs{6.0, 3.0, 0.0};

    nova::linalg::GmresOptions options;
    options.restart = 1;
    options.solver.max_iterations = 1;
    options.solver.relative_tolerance = 1.0e-14;
    options.solver.absolute_tolerance = 1.0e-16;

    const auto result = nova::linalg::solve_gmres(matrix, rhs, options);

    NOVA_REQUIRE(!result.converged);
    NOVA_REQUIRE(result.iterations == 1);
    NOVA_REQUIRE(result.message.find("maximum") != std::string::npos);
}
