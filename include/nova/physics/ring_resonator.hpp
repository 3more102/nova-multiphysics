#pragma once

#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace nova::optics {

struct RingResonatorResult {
    std::complex<double> through_amplitude;
    std::complex<double> circulating_amplitude;
    double through_power;
    double circulating_power;
    double dissipated_power;
    double power_balance_defect;
};

// Single-bus all-pass ring/racetrack, one traveling mode, lossless directional
// coupler. t = real field self-coupling, a = round-trip field attenuation.
// Both are in [0,1]. phase_rad is the real round-trip propagation phase.
// Fields are normalized to unit incident amplitude; circulation is evaluated
// just after the coupler. No backscattering, gain or nonlinear self-consistency.
[[nodiscard]] inline RingResonatorResult all_pass_ring(
    double t, double a, double phase_rad) {
    if (!std::isfinite(t) || !std::isfinite(a) || !std::isfinite(phase_rad) ||
        t < 0.0 || t > 1.0 || a < 0.0 || a > 1.0) {
        throw std::invalid_argument("invalid passive ring parameters");
    }
    // An exactly disconnected ring has no driven circulating field, including
    // the degenerate lossless resonance. Initial stored energy is not modeled.
    if (t == 1.0) return {{1.0, 0.0}, {}, 1.0, 0.0, 0.0, 0.0};

    const double p = std::remainder(phase_rad, 2.0 * std::numbers::pi_v<double>);
    const double s = std::sin(0.5 * p);
    // Stable 1 - exp(i*p) even near resonance.
    const std::complex<double> delta{2.0 * s * s, -std::sin(p)};
    const std::complex<double> den = ((1.0 - t) + t * (1.0 - a)) + t * a * delta;
    const auto h = ((t - a) + a * delta) / den;
    const auto b = std::complex<double>{0.0, std::sqrt((1.0 - t) * (1.0 + t))} / den;
    const double through = std::norm(h);
    const double circulating = std::norm(b);
    const double loss = (1.0 - a) * (1.0 + a) * circulating;
    return {h, b, through, circulating, loss, 1.0 - through - loss};
}

} // namespace nova::optics
