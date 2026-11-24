#include "test_framework.hpp"

#include "nova/linalg/csr_matrix.hpp"

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

NOVA_TEST(csr_accepts_only_canonical_row_storage) {
    using nova::linalg::CsrMatrix;
    const auto m = CsrMatrix::from_sorted_rows(3, 4, {0,2,2,3}, {0,3,2}, {1.0,-2.0,4.0});
    NOVA_REQUIRE(m.exactly_equal(CsrMatrix::from_triplets(3,4,{{0,0,1},{0,3,-2},{2,2,4}})));
    const auto invalid = [](std::vector<std::size_t> offsets, std::vector<std::size_t> columns,
                            std::vector<double> values) {
        try { (void)CsrMatrix::from_sorted_rows(2,2,std::move(offsets),std::move(columns),std::move(values)); }
        catch (const std::invalid_argument&) { return true; }
        return false;
    };
    NOVA_REQUIRE(invalid({0,2,1},{0},{1}));
    NOVA_REQUIRE(invalid({1,1,1},{0},{1}));
    NOVA_REQUIRE(invalid({0,0},{},{}));
    NOVA_REQUIRE(invalid({0,0,2},{0},{1}));
    NOVA_REQUIRE(invalid({0,1,1},{},{1}));
    NOVA_REQUIRE(invalid({0,2,2},{1,0},{1,2}));
    NOVA_REQUIRE(invalid({0,2,2},{1,1},{1,2}));
    NOVA_REQUIRE(invalid({0,1,1},{2},{1}));
    NOVA_REQUIRE(invalid({0,1,1},{0},{0}));
    NOVA_REQUIRE(invalid({0,1,1},{0},{std::numeric_limits<double>::infinity()}));
}
