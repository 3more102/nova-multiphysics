#pragma once
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace nova::optics {
// Lossless continuous-wave plane waves, SI units, no incident second harmonic.
struct ShgInput {
    double vacuum_wavelength_m;
    double fundamental_index;
    double harmonic_index;
    double effective_d_m_per_v;
    double pump_intensity_w_per_m2;
    double length_m;
};

struct ShgResult {
    double harmonic_intensity_w_per_m2;
    double conversion_fraction;
};

namespace shg_detail {
inline long double strength(const ShgInput& p) {
    if (!(std::isfinite(p.vacuum_wavelength_m) && p.vacuum_wavelength_m > 0 &&
          std::isfinite(p.fundamental_index) && p.fundamental_index > 0 &&
          std::isfinite(p.harmonic_index) && p.harmonic_index > 0 &&
          std::isfinite(p.effective_d_m_per_v) &&
          std::isfinite(p.pump_intensity_w_per_m2) && p.pump_intensity_w_per_m2 >= 0 &&
          std::isfinite(p.length_m) && p.length_m >= 0))
        throw std::invalid_argument("invalid SHG input");
    if (p.length_m == 0 || p.pump_intensity_w_per_m2 == 0 || p.effective_d_m_per_v == 0)
        return 0;
    constexpr long double c = 299792458.L;
    constexpr long double eps0 = 8.8541878128e-12L;
    const auto g = (2 * std::numbers::pi_v<long double> /
        static_cast<long double>(p.vacuum_wavelength_m)) *
        std::sqrt(2.L / (eps0 * c)) *
        (std::abs(static_cast<long double>(p.effective_d_m_per_v)) /
         p.fundamental_index) *
        std::sqrt(static_cast<long double>(p.pump_intensity_w_per_m2) /
                  p.harmonic_index) * p.length_m;
    if (!std::isfinite(g)) throw std::overflow_error("SHG strength overflow");
    return g;
}
inline double narrow(long double x) {
    if (!std::isfinite(x) || x > std::numeric_limits<double>::max() || x < 0)
        throw std::overflow_error("SHG result overflow");
    return static_cast<double>(x);
}
} // namespace shg_detail

// Undepleted-pump approximation. Effective mismatch includes any QPM correction.
// Conversion above unity is intentionally not clipped: it signals model breakdown.
[[nodiscard]] inline ShgResult undepleted_shg(const ShgInput& p, double mismatch_per_m) {
    if (!std::isfinite(mismatch_per_m))
        throw std::invalid_argument("invalid SHG phase mismatch");
    const auto g = shg_detail::strength(p);
    if (g == 0) return {0, 0};
    const long double x = static_cast<long double>(mismatch_per_m) * p.length_m / 2;
    if (!std::isfinite(x)) throw std::overflow_error("SHG phase overflow");
    const auto sinc = x == 0 ? 1.L : std::sin(x) / x;
    const auto amplitude = g * sinc;
    const auto f = amplitude * amplitude;
    return {shg_detail::narrow(f * p.pump_intensity_w_per_m2), shg_detail::narrow(f)};
}

struct DepletedShgResult {
    double pump_intensity_w_per_m2;
    double harmonic_intensity_w_per_m2;
    double conversion_fraction;
};

// Exact effective zero-mismatch solution with unseeded harmonic, no absorption.
[[nodiscard]] inline DepletedShgResult phase_matched_depleted_shg(const ShgInput& p) {
    const auto g = shg_detail::strength(p);
    const auto t = std::tanh(g);
    const auto f = t * t;
    // Stable sech^2, retaining the residual pump even when tanh rounds to one.
    const auto e = std::exp(-2 * g);
    const auto remaining = 4 * e / ((1 + e) * (1 + e));
    return {shg_detail::narrow(remaining * p.pump_intensity_w_per_m2),
            shg_detail::narrow(f * p.pump_intensity_w_per_m2),
            shg_detail::narrow(f)};
}
} // namespace nova::optics
