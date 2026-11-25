#include "test_framework.hpp"

#include "nova/linalg/csr_matrix.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

NOVA_TEST(csr_matrix_combines_duplicate_entries) {
    const std::vector<nova::linalg::Triplet> entries{
        {0, 0, 2.0},
        {0, 0, 3.0},
        {0, 1, -1.0},
        {1, 0, 4.0},
        {1, 1, 0.0}
    };

    const auto matrix = nova::linalg::CsrMatrix::from_triplets(2, 2, entries);

    NOVA_REQUIRE(matrix.rows() == 2);
    NOVA_REQUIRE(matrix.columns() == 2);
    NOVA_REQUIRE(matrix.non_zeros() == 3);
    NOVA_REQUIRE_NEAR(matrix.at(0, 0), 5.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(matrix.at(0, 1), -1.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(matrix.at(1, 0), 4.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(matrix.at(1, 1), 0.0, 1.0e-14);
}

NOVA_TEST(csr_matrix_vector_product_is_correct) {
    const std::vector<nova::linalg::Triplet> entries{
        {0, 0, 4.0}, {0, 1, 1.0},
        {1, 0, 1.0}, {1, 1, 3.0}
    };
    const auto matrix = nova::linalg::CsrMatrix::from_triplets(2, 2, entries);
    const std::vector<double> vector{1.0, 2.0};
    const auto result = matrix.multiply(vector);

    NOVA_REQUIRE_NEAR(result[0], 6.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(result[1], 7.0, 1.0e-14);
}

NOVA_TEST(csr_matrix_rejects_out_of_range_triplet) {
    bool threw = false;
    try {
        const std::vector<nova::linalg::Triplet> entries{{2, 0, 1.0}};
        static_cast<void>(nova::linalg::CsrMatrix::from_triplets(2, 2, entries));
    } catch (const std::out_of_range&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}


NOVA_TEST(csr_linear_combination_merges_sparse_entries) {
    using nova::linalg::CsrMatrix;
    using nova::linalg::Triplet;
    const auto a = CsrMatrix::from_triplets(3, 3, {
        {0, 0, 1.0}, {0, 2, 2.0}, {2, 2, 3.0}
    });
    const auto b = CsrMatrix::from_triplets(3, 3, {
        {0, 1, 4.0}, {0, 2, -1.0}, {1, 1, 5.0}
    });
    const auto c = CsrMatrix::linear_combination(a, 2.0, b, 3.0);
    NOVA_REQUIRE_NEAR(c.at(0, 0), 2.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(c.at(0, 1), 12.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(c.at(0, 2), 1.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(c.at(1, 1), 15.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(c.at(2, 2), 6.0, 1.0e-14);
    std::size_t visited = 0;
    c.for_each_nonzero([&](std::size_t, std::size_t, double) {
        ++visited;
    });
    NOVA_REQUIRE(visited == c.non_zeros());
}


NOVA_TEST(csr_exact_equality_checks_dimensions_pattern_and_coefficients) {
    using nova::linalg::CsrMatrix;
    const auto a = CsrMatrix::from_triplets(3, 3, {
        {0,0,4}, {0,1,-1}, {1,0,-1}, {1,1,5}, {2,2,6}
    });
    const auto reorder = CsrMatrix::from_triplets(3, 3, {
        {2,2,6}, {1,1,5}, {1,0,-1}, {0,1,-1}, {0,0,4}
    });
    const auto changed_value = CsrMatrix::from_triplets(3, 3, {
        {0,0,4}, {0,1,-1}, {1,0,-1}, {1,1,5.000000000000001}, {2,2,6}
    });
    const auto changed_pattern = CsrMatrix::from_triplets(3, 3, {
        {0,0,4}, {0,2,-1}, {1,0,-1}, {1,1,5}, {2,2,6}
    });
    const auto changed_dimensions = CsrMatrix::from_triplets(4, 4, {
        {0,0,4}, {0,1,-1}, {1,0,-1}, {1,1,5}, {2,2,6}
    });
    NOVA_REQUIRE(a.exactly_equal(a));
    NOVA_REQUIRE(a.exactly_equal(reorder));
    NOVA_REQUIRE(reorder.exactly_equal(a));
    NOVA_REQUIRE(!a.exactly_equal(changed_value));
    NOVA_REQUIRE(!a.exactly_equal(changed_pattern));
    NOVA_REQUIRE(!a.exactly_equal(changed_dimensions));
}


NOVA_TEST(csr_matrix_vector_product_recovers_overflowing_products_that_cancel) {
    using nova::linalg::CsrMatrix;
    const double coefficient = std::ldexp(1.0, 800);
    const double vector_scale = std::ldexp(1.0, 400);
    const auto matrix = CsrMatrix::from_triplets(1, 2, {
        {0, 0, coefficient}, {0, 1, -coefficient}
    });
    const auto result = matrix.multiply(std::vector<double>{vector_scale, vector_scale});
    NOVA_REQUIRE(std::isfinite(result[0]));
    NOVA_REQUIRE(result[0] == 0.0);
}

NOVA_TEST(csr_matrix_vector_product_preserves_small_residual_after_overflow_cancellation) {
    using nova::linalg::CsrMatrix;
    const double coefficient = std::ldexp(1.0, 800);
    const double vector_scale = std::ldexp(1.0, 400);
    const auto matrix = CsrMatrix::from_triplets(1, 3, {
        {0, 0, coefficient}, {0, 1, -coefficient}, {0, 2, 1.0}
    });
    const auto result = matrix.multiply(
        std::vector<double>{vector_scale, vector_scale, 1.0});
    NOVA_REQUIRE(result[0] == 1.0);
}

NOVA_TEST(csr_matrix_vector_product_recovers_small_term_after_exact_cancellation) {
    using nova::linalg::CsrMatrix;
    const double large = std::ldexp(1.0, 500);
    const auto matrix = CsrMatrix::from_triplets(1, 3, {
        {0, 0, large}, {0, 1, 1.0}, {0, 2, -large}
    });
    const auto result = matrix.multiply(std::vector<double>{1.0, 1.0, 1.0});
    NOVA_REQUIRE(result[0] == 1.0);
}

NOVA_TEST(csr_matrix_vector_product_preserves_representable_subnormal_result) {
    using nova::linalg::CsrMatrix;
    const auto matrix = CsrMatrix::from_triplets(1, 1, {
        {0, 0, std::ldexp(1.0, -600)}
    });
    const auto result = matrix.multiply(std::vector<double>{std::ldexp(1.0, -500)});
    NOVA_REQUIRE(result[0] == std::ldexp(1.0, -1100));
}

NOVA_TEST(csr_matrix_vector_product_keeps_ieee_propagation_for_nonfinite_input) {
    using nova::linalg::CsrMatrix;
    const auto matrix = CsrMatrix::from_triplets(1, 2, {
        {0, 0, 1.0}, {0, 1, 1.0}
    });
    const double infinity = std::numeric_limits<double>::infinity();
    const auto result = matrix.multiply(std::vector<double>{infinity, -infinity});
    NOVA_REQUIRE(std::isnan(result[0]));
}


NOVA_TEST(csr_matrix_vector_product_returns_infinity_for_unrepresentable_exact_result) {
    using nova::linalg::CsrMatrix;
    const auto matrix = CsrMatrix::from_triplets(1, 1, {
        {0, 0, std::ldexp(1.0, 800)}
    });
    const auto result = matrix.multiply(std::vector<double>{std::ldexp(1.0, 400)});
    NOVA_REQUIRE(std::isinf(result[0]));
    NOVA_REQUIRE(result[0] > 0.0);
}
