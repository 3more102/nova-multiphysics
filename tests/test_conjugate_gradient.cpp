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
