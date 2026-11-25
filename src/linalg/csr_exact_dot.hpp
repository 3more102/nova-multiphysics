#pragma once

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>

namespace nova::linalg::detail {

static_assert(sizeof(double) == sizeof(std::uint64_t));
static_assert(std::numeric_limits<double>::is_iec559);
static_assert(std::numeric_limits<double>::digits == 53);

// A finite binary64 product is an integer significand (at most 106 bits)
// multiplied by a power of two. The product exponent range is [-2148, 1942].
// 134 32-bit limbs also accommodate the maximum size_t number of terms.
constexpr int kProductBaseExponent = -2148;
constexpr std::size_t kAccumulatorLimbs = 134;
using ExactAccumulator = std::array<std::uint32_t, kAccumulatorLimbs>;

struct BinaryOperand {
    std::uint64_t significand{};
    int exponent{};
    bool negative{};
};

[[nodiscard]] inline BinaryOperand decode_finite_binary64(const double value) noexcept {
    const std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
    constexpr std::uint64_t fraction_mask = (std::uint64_t{1} << 52U) - 1U;
    const auto exponent_field = static_cast<unsigned>((bits >> 52U) & 0x7ffU);
    const std::uint64_t fraction = bits & fraction_mask;
    if (exponent_field == 0U) {
        return {fraction, -1074, (bits >> 63U) != 0U};
    }
    return {(std::uint64_t{1} << 52U) | fraction,
            static_cast<int>(exponent_field) - 1075,
            (bits >> 63U) != 0U};
}

[[nodiscard]] inline std::array<std::uint32_t, 4> multiply_significands(
    const std::uint64_t left, const std::uint64_t right
) noexcept {
    const std::array<std::uint32_t, 2> a{
        static_cast<std::uint32_t>(left), static_cast<std::uint32_t>(left >> 32U)};
    const std::array<std::uint32_t, 2> b{
        static_cast<std::uint32_t>(right), static_cast<std::uint32_t>(right >> 32U)};
    std::array<std::uint32_t, 4> product{};
    for (std::size_t i = 0; i < 2; ++i) {
        std::uint64_t carry = 0;
        for (std::size_t j = 0; j < 2; ++j) {
            const std::uint64_t current =
                static_cast<std::uint64_t>(a[i]) * b[j] + product[i + j] + carry;
            product[i + j] = static_cast<std::uint32_t>(current);
            carry = current >> 32U;
        }
        std::size_t index = i + 2;
        while (carry != 0U && index < product.size()) {
            const std::uint64_t current =
                static_cast<std::uint64_t>(product[index]) + carry;
            product[index] = static_cast<std::uint32_t>(current);
            carry = current >> 32U;
            ++index;
        }
    }
    return product;
}

inline void add_word(ExactAccumulator& accumulator, std::size_t index,
                     const std::uint32_t word) {
    std::uint64_t carry = word;
    while (carry != 0U && index < accumulator.size()) {
        const std::uint64_t current =
            static_cast<std::uint64_t>(accumulator[index]) + carry;
        accumulator[index] = static_cast<std::uint32_t>(current);
        carry = current >> 32U;
        ++index;
    }
    if (carry != 0U) {
        throw std::overflow_error("CSR exact dot-product accumulator capacity exceeded.");
    }
}

inline void add_shifted_product(ExactAccumulator& accumulator,
                                const std::array<std::uint32_t, 4>& product,
                                const std::size_t bit_shift) {
    const std::size_t word_shift = bit_shift / 32U;
    const unsigned intra_word_shift = static_cast<unsigned>(bit_shift % 32U);
    std::array<std::uint32_t, 5> shifted{};
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < product.size(); ++i) {
        const std::uint64_t current =
            (static_cast<std::uint64_t>(product[i]) << intra_word_shift) | carry;
        shifted[i] = static_cast<std::uint32_t>(current);
        carry = current >> 32U;
    }
    shifted[4] = static_cast<std::uint32_t>(carry);
    for (std::size_t i = 0; i < shifted.size(); ++i) {
        add_word(accumulator, word_shift + i, shifted[i]);
    }
}

inline void add_exact_product(ExactAccumulator& positive, ExactAccumulator& negative,
                              const double left, const double right) {
    const BinaryOperand a = decode_finite_binary64(left);
    const BinaryOperand b = decode_finite_binary64(right);
    if (a.significand == 0U || b.significand == 0U) return;
    const auto product = multiply_significands(a.significand, b.significand);
    const int exponent = a.exponent + b.exponent;
    const std::size_t shift = static_cast<std::size_t>(exponent - kProductBaseExponent);
    ExactAccumulator& destination = (a.negative != b.negative) ? negative : positive;
    add_shifted_product(destination, product, shift);
}

[[nodiscard]] inline int compare_accumulators(const ExactAccumulator& left,
                                              const ExactAccumulator& right) noexcept {
    for (std::size_t i = left.size(); i-- > 0U;) {
        if (left[i] < right[i]) return -1;
        if (left[i] > right[i]) return 1;
    }
    return 0;
}

[[nodiscard]] inline ExactAccumulator subtract_accumulators(
    const ExactAccumulator& larger, const ExactAccumulator& smaller
) noexcept {
    ExactAccumulator result{};
    std::uint64_t borrow = 0;
    for (std::size_t i = 0; i < result.size(); ++i) {
        const std::uint64_t lhs = larger[i];
        const std::uint64_t rhs = static_cast<std::uint64_t>(smaller[i]) + borrow;
        if (lhs >= rhs) {
            result[i] = static_cast<std::uint32_t>(lhs - rhs);
            borrow = 0;
        } else {
            result[i] = static_cast<std::uint32_t>((std::uint64_t{1} << 32U) + lhs - rhs);
            borrow = 1;
        }
    }
    return result;
}

[[nodiscard]] inline bool get_bit(const ExactAccumulator& value,
                                  const std::size_t bit) noexcept {
    const std::size_t word = bit / 32U;
    if (word >= value.size()) return false;
    return ((value[word] >> (bit % 32U)) & 1U) != 0U;
}

[[nodiscard]] inline bool any_bits_below(const ExactAccumulator& value,
                                         const std::size_t bit_count) noexcept {
    const std::size_t full_words = std::min(bit_count / 32U, value.size());
    for (std::size_t i = 0; i < full_words; ++i) {
        if (value[i] != 0U) return true;
    }
    const unsigned remaining = static_cast<unsigned>(bit_count % 32U);
    if (remaining != 0U && full_words < value.size()) {
        const std::uint32_t mask = (std::uint32_t{1} << remaining) - 1U;
        if ((value[full_words] & mask) != 0U) return true;
    }
    return false;
}

[[nodiscard]] inline std::uint64_t extract_bits(const ExactAccumulator& value,
                                                const std::size_t first_bit,
                                                const unsigned count) noexcept {
    std::uint64_t result = 0;
    for (unsigned i = 0; i < count; ++i) {
        if (get_bit(value, first_bit + i)) result |= std::uint64_t{1} << i;
    }
    return result;
}

[[nodiscard]] inline int highest_set_bit(const ExactAccumulator& value) noexcept {
    for (std::size_t i = value.size(); i-- > 0U;) {
        if (value[i] != 0U) {
            return static_cast<int>(i * 32U + 31U - std::countl_zero(value[i]));
        }
    }
    return -1;
}

[[nodiscard]] inline double round_exact_accumulator(const ExactAccumulator& magnitude,
                                                    const bool negative) noexcept {
    const int highest = highest_set_bit(magnitude);
    if (highest < 0) return 0.0;
    const int exact_exponent = highest + kProductBaseExponent;
    const bool subnormal = exact_exponent < -1022;
    const std::size_t shift = subnormal
        ? static_cast<std::size_t>(-1074 - kProductBaseExponent)
        : static_cast<std::size_t>(highest - 52);
    std::uint64_t significand = extract_bits(magnitude, shift, 53U);
    const bool guard = shift > 0U && get_bit(magnitude, shift - 1U);
    const bool sticky = shift > 1U && any_bits_below(magnitude, shift - 1U);
    if (guard && (sticky || (significand & 1U) != 0U)) ++significand;
    if (significand == 0U) return negative ? -0.0 : 0.0;
    const int output_exponent = subnormal ? -1074 :
        kProductBaseExponent + static_cast<int>(shift);
    const double result = std::scalbn(static_cast<double>(significand), output_exponent);
    return negative ? -result : result;
}

// Exact finite-input dot product, rounded once to binary64 (ties to even).
[[nodiscard]] inline double exact_finite_row_dot(
    const std::span<const double> values,
    const std::span<const std::size_t> columns,
    const std::span<const double> vector
) {
    ExactAccumulator positive{};
    ExactAccumulator negative{};
    for (std::size_t i = 0; i < values.size(); ++i) {
        add_exact_product(positive, negative, values[i], vector[columns[i]]);
    }
    const int comparison = compare_accumulators(positive, negative);
    if (comparison == 0) return 0.0;
    const bool is_negative = comparison < 0;
    const ExactAccumulator magnitude = comparison > 0
        ? subtract_accumulators(positive, negative)
        : subtract_accumulators(negative, positive);
    return round_exact_accumulator(magnitude, is_negative);
}

} // namespace nova::linalg::detail
