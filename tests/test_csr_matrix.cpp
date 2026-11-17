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
