#pragma once
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace nova::optics {
struct RingResponse {
    std::complex<double> through;
    std::complex<double> drop;
    double through_power;
    double drop_power;
    double absorbed_power;
};

// Symmetric half-round-trip propagation; lossless real directional couplers.
// t1, t2: self-coupling field amplitudes in [0,1]; a: full-round-trip
// propagation field attenuation in [0,1]. Phase is in radians.
// t2=1 gives an all-pass ring. Unit input at bus 1; no input at bus 2.
// exp(-i omega t) convention, cross-coupling amplitude +i sqrt(1-t^2).
[[nodiscard]] inline RingResponse microring_response(double t1, double t2,
                                                     double a, double phase) {
    const auto valid = [](double x) {
        return std::isfinite(x) && x >= 0.0 && x <= 1.0;
    };
    if (!valid(t1) || !valid(t2) || !valid(a) || !std::isfinite(phase))
        throw std::invalid_argument("invalid passive microring parameters");
    // An uncoupled input cannot excite the ring, even at a singular
    // lossless resonance. Return the externally observable response.
    if (t1 == 1.0) return {{1.0, 0.0}, {}, 1.0, 0.0, 0.0};

    const double half = 0.5 * std::remainder(phase, 4.0 * std::numbers::pi_v<double>);
    const double s = std::sin(half), c = std::cos(half);
    const double b = t2 * a;
    const double q = t1 * b;
    // Difference-of-products form avoids cancellation at high finesse.
    const double one_minus_b = (1.0 - t2) + t2 * (1.0 - a);
    const double one_minus_q = (1.0 - t1) + t1 * one_minus_b;
    const std::complex<double> den{one_minus_q + 2.0*q*s*s, -2.0*q*s*c};
    const double k1 = std::sqrt((1.0-t1)*(1.0+t1));
    const double k2 = std::sqrt((1.0-t2)*(1.0+t2));
    const std::complex<double> numerator{
        (t1-t2) + t2*(1.0-a) + 2.0*b*s*s, -2.0*b*s*c};
    const auto through = numerator / den;
    const auto drop = (-k1*k2*std::sqrt(a)*std::complex<double>{c,s}) / den;
    // Power dissipated over the two half trips, rather than 1-T-D,
    // which loses precision for very small intrinsic absorption.
    const double circulating = std::norm(std::complex<double>{k1,0.0}/den);
    const double absorbed = (1.0-a)*(1.0+t2*t2*a)*circulating;
    return {through, drop, std::norm(through), std::norm(drop), absorbed};
}
} // namespace nova::optics
