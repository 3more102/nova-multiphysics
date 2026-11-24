#pragma once
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace nova::optics {
// Ideal reciprocal, lossless directional couplers and equal half-ring losses.
// t1,t2: real field self-coupling in [0,1]; a: round-trip field survival.
// phase: round-trip propagation phase in radians. Only the input bus is driven.
struct RingResponse {
    std::complex<double> through, drop;
    double through_power, drop_power, absorbed_power;
    double circulating_power; // just after input coupler / incident power
};

[[nodiscard]] inline RingResponse add_drop_ring(double t1, double t2,
                                                double a, double phase) {
    for (double x : {t1,t2,a})
        if (!std::isfinite(x) || x < 0.0 || x > 1.0)
            throw std::invalid_argument("ring coefficients must be in [0,1]");
    if (!std::isfinite(phase))
        throw std::invalid_argument("ring phase must be finite");
    if (t1 == 1.0) return {1.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    const double p = std::remainder(phase, 4.0 * std::numbers::pi_v<double>);
    const double s = std::sin(p / 2.0);
    const double one_minus_cos = 2.0 * s * s;
    const double q = a*t2;
    const double r = t1*q;
    // Preserve tiny coupling/loss terms rather than subtracting a rounded product.
    const double one_minus_q = (1.0-a) + a*(1.0-t2);
    const double one_minus_r = (1.0-t1) + t1*one_minus_q;
    const std::complex<double> den{one_minus_r+r*one_minus_cos, -r*std::sin(p)};
    const double k1sq = (1.0-t1)*(1.0+t1);
    const double k2sq = (1.0-t2)*(1.0+t2);
    const auto h = std::complex<double>{t1-q+q*one_minus_cos, -q*std::sin(p)}/den;
    const auto d = -std::sqrt(k1sq*k2sq*a)*std::polar(1.0,p/2.0)/den;
    const double build = k1sq/std::norm(den);
    const double absorbed = (1.0-a)*(1.0+a*t2*t2)*build;
    const double T = std::norm(h), D = std::norm(d);
    if (!std::isfinite(build) || !std::isfinite(T) || !std::isfinite(D) ||
        !std::isfinite(absorbed) || std::abs(T+D+absorbed-1.0)>1e-8)
        throw std::runtime_error("ring power balance or numerical range failure");
    return {h,d,T,D,absorbed,build};
}
} // namespace nova::optics
