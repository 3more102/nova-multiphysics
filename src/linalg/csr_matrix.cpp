#include "nova/linalg/csr_matrix.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace nova::linalg {

CsrMatrix CsrMatrix::from_triplets(
    const std::size_t rows,
    const std::size_t columns,
    const std::vector<Triplet>& entries,
    const double drop_tolerance
) {
    if (!std::isfinite(drop_tolerance) || drop_tolerance < 0.0) {
        throw std::invalid_argument("CSR drop tolerance must be finite and non-negative.");
    }

    std::vector<std::map<std::size_t, double>> row_entries(rows);

    for (const Triplet& entry : entries) {
        if (entry.row >= rows || entry.column >= columns) {
            throw std::out_of_range("CSR triplet index is outside matrix dimensions.");
        }
        if (!std::isfinite(entry.value)) {
            throw std::invalid_argument("CSR triplet value must be finite.");
        }

        double& accumulated = row_entries[entry.row][entry.column];
        accumulated += entry.value;
        if (!std::isfinite(accumulated)) {
            throw std::overflow_error("CSR duplicate accumulation produced a non-finite value.");
        }
    }

    CsrMatrix matrix;
    matrix.rows_ = rows;
    matrix.columns_ = columns;
    matrix.row_offsets_.resize(rows + 1, 0);

    for (std::size_t row = 0; row < rows; ++row) {
        for (const auto& [column, value] : row_entries[row]) {
            if (std::abs(value) > drop_tolerance) {
                matrix.column_indices_.push_back(column);
                matrix.values_.push_back(value);
            }
        }
        matrix.row_offsets_[row + 1] = matrix.values_.size();
    }

    return matrix;
}

std::size_t CsrMatrix::rows() const noexcept {
    return rows_;
}

std::size_t CsrMatrix::columns() const noexcept {
    return columns_;
}

std::size_t CsrMatrix::non_zeros() const noexcept {
    return values_.size();
}

bool CsrMatrix::exactly_equal(const CsrMatrix& other) const noexcept {
    return rows_ == other.rows_ && columns_ == other.columns_ &&
        row_offsets_ == other.row_offsets_ &&
        column_indices_ == other.column_indices_ && values_ == other.values_;
}

std::vector<double> CsrMatrix::multiply(const std::span<const double> vector) const {
    if (vector.size() != columns_) {
        throw std::invalid_argument("CSR matrix-vector dimension mismatch.");
    }

    std::vector<double> result(rows_, 0.0);
    for (std::size_t row = 0; row < rows_; ++row) {
        double sum = 0.0;
        for (std::size_t index = row_offsets_[row]; index < row_offsets_[row + 1]; ++index) {
            sum += values_[index] * vector[column_indices_[index]];
        }
        result[row] = sum;
    }
    return result;
}

double CsrMatrix::at(const std::size_t row, const std::size_t column) const {
    if (row >= rows_ || column >= columns_) {
        throw std::out_of_range("CSR matrix index is outside matrix dimensions.");
    }

    const auto begin = column_indices_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row]);
    const auto end = column_indices_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row + 1]);
    const auto iterator = std::lower_bound(begin, end, column);

    if (iterator == end || *iterator != column) {
        return 0.0;
    }

    const auto offset = static_cast<std::size_t>(std::distance(column_indices_.begin(), iterator));
    return values_[offset];
}

std::vector<double> CsrMatrix::diagonal() const {
    const std::size_t diagonal_size = std::min(rows_, columns_);
    std::vector<double> result(diagonal_size, 0.0);

    for (std::size_t index = 0; index < diagonal_size; ++index) {
        result[index] = at(index, index);
    }
    return result;
}

CsrMatrix CsrMatrix::linear_combination(
    const CsrMatrix& a, const double scale_a,
    const CsrMatrix& b, const double scale_b
) {
    if (a.rows_ != b.rows_ || a.columns_ != b.columns_) {
        throw std::invalid_argument("CSR linear combination dimension mismatch.");
    }
    if (!std::isfinite(scale_a) || !std::isfinite(scale_b)) {
        throw std::invalid_argument("CSR linear combination scales must be finite.");
    }
    CsrMatrix result;
    result.rows_ = a.rows_;
    result.columns_ = a.columns_;
    result.row_offsets_.reserve(result.rows_ + 1U);
    result.row_offsets_.push_back(0U);
    for (std::size_t row = 0; row < result.rows_; ++row) {
        std::size_t ia = a.row_offsets_[row];
        std::size_t ib = b.row_offsets_[row];
        const std::size_t ea = a.row_offsets_[row + 1U];
        const std::size_t eb = b.row_offsets_[row + 1U];
        while (ia < ea || ib < eb) {
            std::size_t column = 0;
            double value = 0.0;
            if (ib == eb || (ia < ea && a.column_indices_[ia] < b.column_indices_[ib])) {
                column = a.column_indices_[ia];
                value = scale_a * a.values_[ia++];
            } else if (ia == ea || b.column_indices_[ib] < a.column_indices_[ia]) {
                column = b.column_indices_[ib];
                value = scale_b * b.values_[ib++];
            } else {
                column = a.column_indices_[ia];
                value = scale_a * a.values_[ia++] + scale_b * b.values_[ib++];
            }
            if (!std::isfinite(value)) {
                throw std::overflow_error("CSR linear combination produced a non-finite value.");
            }
            if (value != 0.0) {
                result.column_indices_.push_back(column);
                result.values_.push_back(value);
            }
        }
        result.row_offsets_.push_back(result.values_.size());
    }
    return result;
}

} // namespace nova::linalg
