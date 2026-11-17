#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace nova::linalg {

struct Triplet {
    std::size_t row{};
    std::size_t column{};
    double value{};
};

class CsrMatrix {
public:
    [[nodiscard]] static CsrMatrix from_triplets(
        std::size_t rows,
        std::size_t columns,
        const std::vector<Triplet>& entries,
        double drop_tolerance = 0.0
    );

    [[nodiscard]] std::size_t rows() const noexcept;
    [[nodiscard]] std::size_t columns() const noexcept;
    [[nodiscard]] std::size_t non_zeros() const noexcept;

    [[nodiscard]] std::vector<double> multiply(
        std::span<const double> vector
    ) const;

    [[nodiscard]] double at(std::size_t row, std::size_t column) const;
    [[nodiscard]] std::vector<double> diagonal() const;

private:
    std::size_t rows_{};
    std::size_t columns_{};
    std::vector<std::size_t> row_offsets_;
    std::vector<std::size_t> column_indices_;
    std::vector<double> values_;
};

} // namespace nova::linalg
