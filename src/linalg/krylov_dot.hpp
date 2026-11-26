#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace nova::linalg::detail {

// A mantissa/exponent inner product prevents tiny Krylov products from
// underflowing when their ratio (alpha or beta) is representable.
struct ScaledDot { double mantissa{}; int exponent{}; };

[[nodiscard]] inline ScaledDot dot(const std::vector<double>& left, const std::vector<double>& right) {
    double direct = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) direct += left[i] * right[i];
    if (std::isnormal(direct)) {
        int exponent = 0;
        return {std::frexp(direct, &exponent), exponent};
    }
    double left_scale = 0.0, right_scale = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        left_scale = std::max(left_scale, std::abs(left[i]));
        right_scale = std::max(right_scale, std::abs(right[i]));
    }
    if (left_scale == 0.0 || right_scale == 0.0) return {};
    double normalized = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        normalized += (left[i] / left_scale) * (right[i] / right_scale);
    }
    if (normalized == 0.0 || !std::isfinite(normalized)) return {normalized, 0};
    int left_exp = 0, right_exp = 0, sum_exp = 0;
    const double lm = std::frexp(left_scale, &left_exp);
    const double rm = std::frexp(right_scale, &right_exp);
    const double sm = std::frexp(normalized, &sum_exp);
    return {lm * rm * sm, left_exp + right_exp + sum_exp};
}

[[nodiscard]] inline double dot_ratio(const ScaledDot& numerator, const ScaledDot& denominator) {
    return std::scalbn(numerator.mantissa / denominator.mantissa,
                       numerator.exponent - denominator.exponent);
}

// Preserve the range of ||u|| ||v|| in a dimensionless breakdown test.
[[nodiscard]] inline ScaledDot norm_product(double left, double right) {
    int left_exponent = 0, right_exponent = 0;
    const double left_mantissa = std::frexp(left, &left_exponent);
    const double right_mantissa = std::frexp(right, &right_exponent);
    return {left_mantissa * right_mantissa, left_exponent + right_exponent};
}

} // namespace nova::linalg::detail
