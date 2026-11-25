#include "test_framework.hpp"

#include "nova/linalg/bicgstab.hpp"

#include <cmath>
#include <limits>
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

NOVA_TEST(bicgstab_preserves_nonsymmetric_solution_across_binary_scales) {
    using namespace nova::linalg;
    SolverOptions options;
    options.absolute_tolerance = 0.0;
    options.relative_tolerance = 1.0e-12;

    // [[4,1],[-2,3]] (1,2)^T = (6,4)^T. Vary the load alone and
    // the entire equation; in each case start from a nonzero physical guess.
    for (int exponent : {-1000, -660, -540, 0, 60, 540, 660, 1000}) {
        const double load_scale = std::scalbn(1.0, exponent);
        for (double matrix_scale : {1.0, load_scale}) {
            const double field_scale = load_scale / matrix_scale;
            const auto matrix = CsrMatrix::from_triplets(2, 2, {
                {0, 0, 4.0 * matrix_scale}, {0, 1, matrix_scale},
                {1, 0, -2.0 * matrix_scale}, {1, 1, 3.0 * matrix_scale}
            });
            const std::vector<double> rhs{6.0 * load_scale, 4.0 * load_scale};
            const std::vector<double> initial{0.5 * field_scale, -0.5 * field_scale};
            for (auto method : {PreconditionerType::None, PreconditionerType::Jacobi}) {
                options.preconditioner = method;
                const auto result = solve_bicgstab(matrix, rhs, options, initial);
                NOVA_REQUIRE(result.converged);
                NOVA_REQUIRE(result.iterations > 0);
                const double x = result.solution[0] / field_scale;
                const double y = result.solution[1] / field_scale;
                NOVA_REQUIRE_NEAR(x, 1.0, 1.0e-12);
                NOVA_REQUIRE_NEAR(y, 2.0, 1.0e-12);
                NOVA_REQUIRE_NEAR(result.initial_residual_norm / load_scale,
                                  std::hypot(4.5, 6.5), 2.0e-15);
                const double residual = std::hypot(6.0 - (4.0 * x + y), 4.0 - (-2.0 * x + 3.0 * y));
                NOVA_REQUIRE_NEAR(result.residual_norm / load_scale, residual, 2.0e-15);
                NOVA_REQUIRE(residual <= options.relative_tolerance * std::hypot(4.5, 6.5));

                const auto repeat = solve_bicgstab(matrix, rhs, options, initial);
                NOVA_REQUIRE(repeat.solution == result.solution);
                NOVA_REQUIRE(repeat.iterations == result.iterations);
                NOVA_REQUIRE(repeat.residual_norm == result.residual_norm);
            }
        }
    }
}

NOVA_TEST(bicgstab_confirms_recurrence_convergence_after_residual_replacement) {
    using namespace nova::linalg;
    const double scale = std::scalbn(1.0, -660);
    const auto matrix = CsrMatrix::from_triplets(2, 2, {
        {0, 0, 4.0 * scale}, {0, 1, scale}, {1, 0, -2.0 * scale}, {1, 1, 3.0 * scale}
    });
    const std::vector<double> rhs{6.0 * scale, 4.0 * scale};
    SolverOptions options;
    options.absolute_tolerance = 0.0;
    options.relative_tolerance = 1.0e-30;
    options.max_iterations = 32;
    // This dyadic load admits a zero residual. Roundoff can make the
    // recursive residual meet the target before the returned field does.
    const auto result = solve_bicgstab(matrix, rhs, options);
    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.iterations > 0);
    const double x = result.solution[0], y = result.solution[1];
    NOVA_REQUIRE(4.0 * x + y == 6.0);
    NOVA_REQUIRE(-2.0 * x + 3.0 * y == 4.0);
    NOVA_REQUIRE(result.residual_norm == 0.0);
}

NOVA_TEST(bicgstab_preserves_genuine_breakdown_across_scales) {
    using namespace nova::linalg;
    SolverOptions options;
    options.absolute_tolerance = 0.0;
    for (int exponent : {-660, 0, 660}) {
        const double scale = std::scalbn(1.0, exponent);
        const std::vector<double> rhs{scale, 0.0};
        // A nonsingular system can still break down: after alpha, s is
        // parallel to e2 while A*s is parallel to e1, so omega is zero.
        const auto matrix = CsrMatrix::from_triplets(2, 2, {
            {0, 0, scale}, {0, 1, 2.0 * scale}, {1, 0, -scale}
        });
        const auto result = solve_bicgstab(matrix, rhs, options);
        NOVA_REQUIRE(!result.converged);
        NOVA_REQUIRE(result.message.find("orthogonal") != std::string::npos);
        NOVA_REQUIRE(result.initial_residual_norm == scale);

        const auto zero = CsrMatrix::from_triplets(2, 2, {});
        const auto singular = solve_bicgstab(zero, rhs, options);
        NOVA_REQUIRE(!singular.converged);
        NOVA_REQUIRE(singular.iterations == 0);
    }
}

NOVA_TEST(bicgstab_rejects_unrepresentable_initial_residual_norm) {
    using namespace nova::linalg;
    const auto matrix = CsrMatrix::from_triplets(2, 2, {{0, 0, 1.0}, {1, 1, 1.0}});
    const std::vector<double> rhs(2, std::numeric_limits<double>::max());
    bool threw = false;
    try {
        static_cast<void>(solve_bicgstab(matrix, rhs));
    } catch (const std::overflow_error&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}
