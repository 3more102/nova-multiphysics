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

    // Take ownership of assembled canonical rows without rebuilding via triplets.
    // Offsets span all entries, columns increase strictly per row, and stored
    // coefficients must be finite and nonzero.
    [[nodiscard]] static CsrMatrix from_sorted_rows(
        std::size_t rows, std::size_t columns, std::vector<std::size_t> row_offsets,
        std::vector<std::size_t> column_indices, std::vector<double> values
    );

    [[nodiscard]] std::size_t rows() const noexcept;
    [[nodiscard]] std::size_t columns() const noexcept;
    [[nodiscard]] std::size_t non_zeros() const noexcept;

    // Exact identity of canonical CSR dimensions, sparsity pattern and values.
    // Compares the storage arrays directly in O(rows + nnz).
    [[nodiscard]] bool exactly_equal(const CsrMatrix& other) const noexcept;

    [[nodiscard]] std::vector<double> multiply(
        std::span<const double> vector
    ) const;

    [[nodiscard]] double at(std::size_t row, std::size_t column) const;
    [[nodiscard]] std::vector<double> diagonal() const;

    // Principal block [first,first+count) with diagonal left/right scaling.
    // Scale arrays are full-system length. Traverses only the selected rows.
    [[nodiscard]] CsrMatrix scaled_block(std::size_t first, std::size_t count,
        std::span<const double> left, std::span<const double> right) const;

    // Sparse linear combination; merges sorted CSR rows in O(nnz(A)+nnz(B)).
    [[nodiscard]] static CsrMatrix linear_combination(
        const CsrMatrix& a, double scale_a,
        const CsrMatrix& b, double scale_b
    );

    template <class Visitor>
    void for_each_nonzero(Visitor&& visit) const {
        for (std::size_t row = 0; row < rows_; ++row) {
            for (std::size_t index = row_offsets_[row];
                 index < row_offsets_[row + 1]; ++index) {
                visit(row, column_indices_[index], values_[index]);
            }
        }
    }

private:
    std::size_t rows_{};
    std::size_t columns_{};
    std::vector<std::size_t> row_offsets_;
    std::vector<std::size_t> column_indices_;
    std::vector<double> values_;
};

// Immutable assembly graph. Numerical zeros can enter and leave the canonical
// matrix on each evaluation without rebuilding or sorting the element graph.
class CsrAssemblyPattern {
public:
    CsrAssemblyPattern() = default;
    CsrAssemblyPattern(std::size_t columns, std::vector<std::size_t> offsets,
        std::vector<std::size_t> indices);
    [[nodiscard]] std::size_t slots() const noexcept { return indices_.size(); }
    [[nodiscard]] std::size_t slot(std::size_t row, std::size_t column) const;
    [[nodiscard]] CsrMatrix matrix(std::span<const double> values) const;
private:
    std::size_t columns_{};
    std::vector<std::size_t> offsets_, indices_;
};

} // namespace nova::linalg
