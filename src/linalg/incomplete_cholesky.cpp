#include "nova/linalg/incomplete_cholesky.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace nova::linalg {

IncompleteCholesky0::IncompleteCholesky0(const CsrMatrix& matrix) {
    const auto n = matrix.rows();
    if (n == 0 || n != matrix.columns()) {
        throw std::invalid_argument("IC(0) requires a non-empty square matrix.");
    }
    constexpr double roundoff = std::numeric_limits<double>::epsilon();
    auto scale = matrix.diagonal();
    for (std::size_t i = 0; i < n; ++i) {
        if (!(scale[i] > 0.0) || !std::isfinite(scale[i])) {
            throw std::invalid_argument("IC(0) requires a positive diagonal at row " + std::to_string(i));
        }
        scale[i] = std::sqrt(scale[i]);
    }

    if (matrix.non_zeros() < n) {
        throw std::invalid_argument("IC(0) matrix is missing diagonal entries.");
    }
    row_offsets_.assign(n + 1, 0);
    const auto off_diagonal_pairs = (matrix.non_zeros() - n) / 2;
    columns_.reserve(off_diagonal_pairs);
    lower_.reserve(off_diagonal_pairs);
    matrix.for_each_nonzero([&](std::size_t i, std::size_t j, double a) {
        if (i == j) return;
        const double transpose = matrix.at(j, i);
        if (!std::isfinite(a) || !std::isfinite(transpose) || transpose == 0.0 || std::abs(a - transpose) >
            64 * roundoff * std::max(std::abs(a), std::abs(transpose))) {
            throw std::invalid_argument("IC(0) requires symmetric CSR entries at row " + std::to_string(i));
        }
        if (i > j) {
            // Diagonally equilibrate A. Dividing by the smaller root first
            // avoids losing a coupling solely to a large disparity of units.
            lower_.push_back((a / std::min(scale[i], scale[j])) / std::max(scale[i], scale[j]));
            columns_.push_back(j);
            ++row_offsets_[i + 1];
        }
    });
    std::partial_sum(row_offsets_.begin(), row_offsets_.end(), row_offsets_.begin());

    diagonal_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        double pivot = 1.0;
        for (auto p = row_offsets_[i]; p < row_offsets_[i+1]; ++p) {
            const auto j = columns_[p];
            double entry = lower_[p];
            // Subtract products only at common, already factored columns.
            auto left = row_offsets_[i], right = row_offsets_[j];
            while (left < p && right < row_offsets_[j+1]) {
                if (columns_[left] < columns_[right]) {
                    ++left;
                } else if (columns_[right] < columns_[left]) {
                    ++right;
                } else {
                    entry = std::fma(-lower_[left], lower_[right], entry);
                    ++left;
                    ++right;
                }
            }
            lower_[p] = entry / diagonal_[j];
            if (!std::isfinite(lower_[p])) {
                throw std::invalid_argument("IC(0) produced a non-finite factor entry.");
            }
            pivot = std::fma(-lower_[p], lower_[p], pivot);
        }
        if (!(pivot > 32 * roundoff) || !std::isfinite(pivot)) {
            throw std::invalid_argument("IC(0) has a non-positive or unresolved pivot at row " + std::to_string(i));
        }
        diagonal_[i] = std::sqrt(pivot);
    }

    // Return the factor to the original variables after all row intersections.
    for (std::size_t i = 0; i < n; ++i) {
        diagonal_[i] *= scale[i];
        if (!std::isfinite(diagonal_[i]) || !(diagonal_[i] > 0.0)) {
            throw std::invalid_argument("IC(0) diagonal factor is outside the finite range.");
        }
        for (auto p = row_offsets_[i]; p < row_offsets_[i+1]; ++p) {
            lower_[p] *= scale[i];
            if (!std::isfinite(lower_[p])) {
                throw std::invalid_argument("IC(0) lower factor is outside the finite range.");
            }
        }
    }
}

std::vector<double> IncompleteCholesky0::apply(std::span<const double> rhs) const {
    if (rhs.size() != dimension() ||
        !std::all_of(rhs.begin(), rhs.end(), [](double x) { return std::isfinite(x); })) {
        throw std::invalid_argument("IC(0) needs a finite right-hand side of matching dimension.");
    }
    std::vector<double> result(rhs.begin(), rhs.end());
    for (std::size_t i = 0; i < dimension(); ++i) {
        double value = result[i];
        for (auto p = row_offsets_[i]; p < row_offsets_[i+1]; ++p) {
            value = std::fma(-lower_[p], result[columns_[p]], value);
        }
        result[i] = value / diagonal_[i];
        if (!std::isfinite(result[i])) throw std::overflow_error("IC(0) forward solve overflow.");
    }
    // Scatter each known component into the preceding rows for L^T.
    for (std::size_t i = dimension(); i-- > 0;) {
        result[i] /= diagonal_[i];
        if (!std::isfinite(result[i])) throw std::overflow_error("IC(0) backward solve overflow.");
        for (auto p = row_offsets_[i]; p < row_offsets_[i+1]; ++p) {
            auto& value = result[columns_[p]];
            value = std::fma(-lower_[p], result[i], value);
            if (!std::isfinite(value)) throw std::overflow_error("IC(0) backward update overflow.");
        }
    }
    return result;
}

} // namespace nova::linalg
