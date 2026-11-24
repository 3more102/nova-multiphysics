#pragma once
#include <cmath>
#include <complex>
#include <stdexcept>

namespace nova::optics {
struct BraggGratingResult {
    std::complex<double> reflection;
    std::complex<double> transmission;
    double reflectance;
    double transmittance;
};

// Uniform lossless two-mode grating: A(0)=1, B(L)=0, r=B(0), t=A(L).
// A'=i*delta*A+i*kappa*B; B'=-i*kappa*A-i*delta*B.
// Real detuning/coupling in m^-1, length in m; envelopes omit carrier phase.
[[nodiscard]] inline BraggGratingResult uniform_bragg_grating(
    double detuning_per_m, double coupling_per_m, double length_m) {
    using C = std::complex<double>;
    if (!std::isfinite(detuning_per_m) || !std::isfinite(coupling_per_m) ||
        !std::isfinite(length_m) || length_m < 0.0)
        throw std::invalid_argument("invalid uniform grating inputs");
    const double d = detuning_per_m * length_m;
    const double k = coupling_per_m * length_m;
    if (!std::isfinite(d) || !std::isfinite(k))
        throw std::overflow_error("nonfinite grating phase");
    if (k == 0.0) return {C{}, std::polar(1.0, d), 0.0, 1.0};
    const double a = std::abs(k), b = std::abs(d);
    const double sum = a + b;
    if (!std::isfinite(sum)) throw std::overflow_error("grating phase sum overflow");
    const double x = std::sqrt(std::abs(a-b)) * std::sqrt(sum);
    if (!std::isfinite(x)) throw std::overflow_error("grating exponent overflow");
    double c, s, scale;
    if (a > b && x > 0.0) {
        // Factor exp(x) out of the transfer matrix; no growing exponential.
        c = 0.5 * (1.0 + std::exp(-2.0*x));
        s = -std::expm1(-2.0*x) / x * 0.5;
        scale = std::exp(-x);
    } else {
        c = std::cos(x);
        s = x == 0.0 ? 1.0 : std::sin(x)/x;
        scale = 1.0;
    }
    const C denominator{c, -d*s};
    const C r = C{0.0, k*s}/denominator;
    const C t = scale/denominator;
    const double R = std::norm(r), T = std::norm(t);
    if (!std::isfinite(R) || !std::isfinite(T) || std::abs(R+T-1.0) > 1e-9)
        throw std::runtime_error("uniform grating power balance failed");
    return {r,t,R,T};
}
} // namespace nova::optics
