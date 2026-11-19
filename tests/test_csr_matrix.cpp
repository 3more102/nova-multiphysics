#include "test_framework.hpp"

#include "nova/linalg/csr_matrix.hpp"

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

NOVA_TEST(csr_exact_entry_comparison_checks_values_pattern_and_dimensions) {
    using nova::linalg::CsrMatrix;
    const auto a = CsrMatrix::from_triplets(3, 4, {
        {2,3,2.0}, {0,2,1.0}, {2,1,-3.0}
    });
    // Input triplet order and cancelled duplicate entries do not affect
    // equality of the canonical stored matrix.
    const auto same = CsrMatrix::from_triplets(3, 4, {
        {2,1,-3.0}, {0,2,0.5}, {2,3,2.0}, {0,2,0.5},
        {1,0,1.0}, {1,0,-1.0}
    });
    const auto new_value = CsrMatrix::from_triplets(3, 4, {
        {0,2,1.0}, {2,1,-3.0}, {2,3,3.0}
    });
    const auto new_column = CsrMatrix::from_triplets(3, 4, {
        {0,2,1.0}, {2,1,-3.0}, {2,2,2.0}
    });
    const auto new_row = CsrMatrix::from_triplets(3, 4, {
        {0,2,1.0}, {1,3,2.0}, {2,1,-3.0}
    });
    const auto new_shape = CsrMatrix::from_triplets(4, 4, {
        {0,2,1.0}, {2,1,-3.0}, {2,3,2.0}
    });
    NOVA_REQUIRE(a.has_identical_entries(same));
    NOVA_REQUIRE(same.has_identical_entries(a));
    NOVA_REQUIRE(!a.has_identical_entries(new_value));
    NOVA_REQUIRE(!a.has_identical_entries(new_column));
    NOVA_REQUIRE(!a.has_identical_entries(new_row));
    NOVA_REQUIRE(!a.has_identical_entries(new_shape));
    NOVA_REQUIRE(CsrMatrix::from_triplets(0,0,{}).has_identical_entries(
        CsrMatrix::from_triplets(0,0,{})));
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
