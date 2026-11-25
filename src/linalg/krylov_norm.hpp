#pragma once

#include <cmath>
#include <span>

namespace nova::linalg::detail {

[[nodiscard]] inline double euclidean_norm(std::span<const double> vector) {
    double squared = 0.0;
    for (double value : vector) squared += value * value;
    if (std::isnormal(squared)) return std::sqrt(squared);
    // Scaled sum of squares (Blue): avoids both squaring underflow and overflow.
    double scale = 0.0, sum = 1.0;
    for (double value : vector) {
        const double magnitude = std::abs(value);
        if (magnitude == 0.0) continue;
        if (magnitude > scale) {
            const double ratio = scale / magnitude;
            sum = 1.0 + sum * ratio * ratio;
            scale = magnitude;
        } else {
            const double ratio = magnitude / scale;
            sum += ratio * ratio;
        }
    }
    return scale * std::sqrt(sum);
}

} // namespace nova::linalg::detail
