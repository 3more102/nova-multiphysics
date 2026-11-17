#include "test_framework.hpp"

#include "nova/nonlinear/newton_solver.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

nova::nonlinear::NonlinearResidualFunction square_root_residual() {
    return [](const std::span<const double> x) {
        return std::vector<double>{x[0] * x[0] - 2.0};
    };
}

nova::nonlinear::NonlinearJacobianFunction square_root_jacobian() {
    return [](const std::span<const double> x) {
        const std::vector<nova::linalg::Triplet> entries{
            {0, 0, 2.0 * x[0]}
        };
        return nova::linalg::CsrMatrix::from_triplets(1, 1, entries);
    };
}

} // namespace

NOVA_TEST(newton_solves_scalar_square_root_problem) {
    const std::vector<double> initial_guess{1.0};

    nova::nonlinear::NewtonOptions options;
    options.residual_tolerance = 1.0e-12;
    options.gmres.solver.relative_tolerance = 1.0e-12;
    options.gmres.solver.absolute_tolerance = 1.0e-14;

    const auto result = nova::nonlinear::solve_newton(
        square_root_residual(),
        square_root_jacobian(),
        initial_guess,
        options
    );

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE_NEAR(result.solution[0], std::sqrt(2.0), 1.0e-10);
    NOVA_REQUIRE(result.residual_norm <= options.residual_tolerance);
    NOVA_REQUIRE(result.iterations > 0);
}

NOVA_TEST(newton_backtracking_rejects_non_finite_trial_and_converges) {
    const std::vector<double> initial_guess{0.1};
    const auto guarded_residual = [](const std::span<const double> x) {
        if (x[0] > 2.0) {
            return std::vector<double>{
                std::numeric_limits<double>::infinity()
            };
        }
        return std::vector<double>{x[0] * x[0] - 2.0};
    };

    nova::nonlinear::NewtonOptions options;
    options.residual_tolerance = 1.0e-10;
    options.gmres.solver.relative_tolerance = 1.0e-12;
    options.gmres.solver.absolute_tolerance = 1.0e-14;

    const auto result = nova::nonlinear::solve_newton(
        guarded_residual,
        square_root_jacobian(),
        initial_guess,
        options
    );

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE_NEAR(result.solution[0], std::sqrt(2.0), 1.0e-8);
}

NOVA_TEST(newton_supports_bicgstab_for_linearized_systems) {
    const std::vector<double> initial_guess{1.0};

    nova::nonlinear::NewtonOptions options;
    options.linear_solver = nova::nonlinear::NewtonLinearSolver::BiCGSTAB;
    options.residual_tolerance = 1.0e-12;
    options.bicgstab.relative_tolerance = 1.0e-12;
    options.bicgstab.absolute_tolerance = 1.0e-14;

    const auto result = nova::nonlinear::solve_newton(
        square_root_residual(),
        square_root_jacobian(),
        initial_guess,
        options
    );

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE_NEAR(result.solution[0], std::sqrt(2.0), 1.0e-10);
}

NOVA_TEST(newton_reports_singular_linearized_system) {
    const std::vector<double> initial_guess{0.0};
    const auto residual = [](const std::span<const double>) {
        return std::vector<double>{1.0};
    };
    const auto jacobian = [](const std::span<const double>) {
        const std::vector<nova::linalg::Triplet> entries{{0, 0, 0.0}};
        return nova::linalg::CsrMatrix::from_triplets(1, 1, entries);
    };

    const auto result = nova::nonlinear::solve_newton(
        residual,
        jacobian,
        initial_guess
    );

    NOVA_REQUIRE(!result.converged);
    NOVA_REQUIRE(result.message.find("linear solve failed") != std::string::npos);
}

NOVA_TEST(newton_rejects_residual_dimension_mismatch) {
    const std::vector<double> initial_guess{1.0, 2.0};
    const auto residual = [](const std::span<const double>) {
        return std::vector<double>{1.0};
    };
    const auto jacobian = [](const std::span<const double>) {
        const std::vector<nova::linalg::Triplet> entries{
            {0, 0, 1.0}, {1, 1, 1.0}
        };
        return nova::linalg::CsrMatrix::from_triplets(2, 2, entries);
    };

    bool threw = false;
    try {
        static_cast<void>(nova::nonlinear::solve_newton(
            residual,
            jacobian,
            initial_guess
        ));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}
