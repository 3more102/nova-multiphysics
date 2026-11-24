#pragma once
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace nova::optics {
struct MicroringResult {
    std::complex<double> through_amplitude, drop_amplitude;
    double through_power, drop_power, absorbed_power;
    double circulating_power, power_balance_defect;
};
// Single input, lossless reciprocal couplers with real field self-transmission.
// Round-trip field attenuation a in [0,1]; equal propagation halves.
// Phase in radians, exp(-i omega t). t2=1 gives an all-pass ring.
[[nodiscard]] inline MicroringResult microring(double t1, double t2, double a,
                                                double phase_rad) {
    const auto unit = [](double x) { return std::isfinite(x) && x>=0 && x<=1; };
    if (!unit(t1) || !unit(t2) || !unit(a) || !std::isfinite(phase_rad))
        throw std::invalid_argument("invalid microring parameters");
    using C = std::complex<double>;
    const double p = std::remainder(phase_rad, 4*std::numbers::pi_v<double>);
    const double s = std::sin(p/2), c = std::cos(p/2);
    const C one_minus_z{2*s*s, -2*s*c};
    const double g = a*t1*t2;
    const C den = (1-g) + g*one_minus_z;
    if (den == C{} || (g == 1 && std::remainder(phase_rad,
        2*std::numbers::pi_v<double>) == 0)) throw std::domain_error("singular lossless uncoupled ring");
    const double k1 = std::sqrt((1-t1)*(1+t1));
    const double k2 = std::sqrt((1-t2)*(1+t2));
    const C through = ((t1-a*t2) + a*t2*one_minus_z)/den;
    const C drop = -k1*k2*std::sqrt(a)*C{c,s}/den;
    const double circulating = std::norm(C{0,k1}/den);
    const double T = std::norm(through), D = std::norm(drop);
    // Independently integrate power lost on both half-round-trip segments.
    const double A = (1-a)*(1+a*t2*t2)*circulating;
    const double defect = 1-T-D-A;
    if (!std::isfinite(circulating) || !std::isfinite(T) || !std::isfinite(D) ||
        !std::isfinite(A) || std::abs(defect)>1e-8)
        throw std::runtime_error("microring power balance failed");
    return {through,drop,T,D,A,circulating,defect};
}
} // namespace nova::optics
