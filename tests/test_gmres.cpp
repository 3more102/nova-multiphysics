#include "test_framework.hpp"

#include "nova/linalg/gmres.hpp"

#include <cmath>
#include <limits>
#include <utility>
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

NOVA_TEST(gmres_preserves_solution_and_true_residual_across_binary_scales) {
    using namespace nova::linalg;
    // Exact solution of the nonsymmetric stencil is amplitude * (1, 2, -1).
    // Change the excitation and the units of the entire equation independently.
    for (const auto& [matrix_exp, solution_exp] : {
             std::pair{0, -660}, {0, 0}, {0, 660}, {-660, 0}, {660, 0}}) {
        const double a = std::ldexp(1.0, matrix_exp);
        const double amplitude = std::ldexp(1.0, solution_exp);
        const auto matrix = CsrMatrix::from_triplets(3, 3, {
            {0, 0, 4*a}, {0, 1, a}, {1, 0, -2*a}, {1, 1, 3*a},
            {1, 2, a}, {2, 1, a}, {2, 2, 2*a}});
        const double load = a * amplitude;
        const std::vector<double> rhs{6*load, 3*load, 0};
        const std::vector<double> guess{0.5*amplitude, 0.0, 0.0};
        for (auto pc : {PreconditionerType::None, PreconditionerType::Jacobi}) {
            GmresOptions options;
            options.restart = 2;
            options.solver.relative_tolerance = 1e-12;
            options.solver.absolute_tolerance = 0;
            options.solver.preconditioner = pc;
            const auto r = solve_gmres(matrix, rhs, options, guess);
            NOVA_REQUIRE(r.converged && r.iterations > options.restart);
            const double x = r.solution[0] / amplitude;
            const double y = r.solution[1] / amplitude;
            const double z = r.solution[2] / amplitude;
            NOVA_REQUIRE_NEAR(x, 1, 4e-12);
            NOVA_REQUIRE_NEAR(y, 2, 4e-12);
            NOVA_REQUIRE_NEAR(z, -1, 4e-12);
            const double initial = std::hypot(4.0, 4.0);
            const double residual = std::hypot(6-(4*x+y), 3-(-2*x+3*y+z), y+2*z);
            NOVA_REQUIRE_NEAR(r.initial_residual_norm / load, initial, 2e-15);
            NOVA_REQUIRE_NEAR(r.residual_norm / load, residual, 3e-15);
            NOVA_REQUIRE(residual <= options.solver.relative_tolerance * initial);
        }
    }
}

NOVA_TEST(gmres_preserves_absolute_tolerance_and_true_norm_overflow) {
    using namespace nova::linalg;
    const auto identity = CsrMatrix::from_triplets(2, 2, {{0,0,1}, {1,1,1}});
    const std::vector<double> tiny{3e-200, 4e-200};
    GmresOptions opts;
    opts.solver.relative_tolerance = 0;
    opts.solver.absolute_tolerance = 6e-200;
    const auto accepted = solve_gmres(identity, tiny, opts);
    NOVA_REQUIRE(accepted.converged && accepted.iterations == 0);
    NOVA_REQUIRE_NEAR(accepted.residual_norm / 1e-200, 5, 2e-15);
    opts.solver.absolute_tolerance = 4e-200;
    const auto solved = solve_gmres(identity, tiny, opts);
    NOVA_REQUIRE(solved.converged && solved.iterations == 1);
    // Each component and its identity solution are exactly representable;
    // the physical norm rounds to one unit of the subnormal grid.
    const double unit = std::numeric_limits<double>::denorm_min();
    opts.solver.absolute_tolerance = 0;
    opts.solver.relative_tolerance = 1e-8;
    const std::vector<double> subnormal{unit, unit};
    const auto exact = solve_gmres(identity, subnormal, opts);
    NOVA_REQUIRE(exact.converged && exact.iterations == 1);
    NOVA_REQUIRE(exact.solution == subnormal && exact.residual_norm == 0);

    const double big = std::numeric_limits<double>::max();
    bool rejected = false;
    try {
        (void)solve_gmres(identity, std::vector<double>{big, big}, opts);
    } catch (const std::overflow_error&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}
