#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <span>
#include <stdexcept>

namespace nova::physics {

// All frequencies and damping rates are angular frequencies in rad/s.
// A zero resonance frequency represents a Drude free-carrier term.
struct OpticalOscillator {
    double plasma_frequency{};
    double resonance_frequency{};
    double damping_rate{};
};

// Passive scalar Drude-Lorentz permittivity, exp(-i omega t) convention.
// epsilon_r = epsilon_inf + sum wp^2/(w0^2-w^2-i gamma w).
// Rejects exact undamped poles and nonrepresentable results.
[[nodiscard]] inline std::complex<double> optical_permittivity(
    double omega, double epsilon_infinity,
    std::span<const OpticalOscillator> oscillators) {
    if (!(omega>0.0) || !std::isfinite(omega)
        || !(epsilon_infinity>0.0) || !std::isfinite(epsilon_infinity))
        throw std::invalid_argument("Invalid optical frequency or background permittivity.");
    std::complex<double> epsilon{epsilon_infinity,0.0};
    for (const auto& o:oscillators) {
        if (!std::isfinite(o.plasma_frequency) || o.plasma_frequency<0.0
            || !std::isfinite(o.resonance_frequency) || o.resonance_frequency<0.0
            || !std::isfinite(o.damping_rate) || o.damping_rate<0.0)
            throw std::invalid_argument("Invalid passive optical oscillator.");
        if (o.plasma_frequency==0.0) continue;
        // Scale before multiplying to avoid squaring optical-scale frequencies.
        const double scale=std::max({omega,o.resonance_frequency,o.damping_rate});
        const double w=omega/scale,r=o.resonance_frequency/scale;
        const double real=((o.resonance_frequency-omega)/scale)*(r+w);
        const double loss=(o.damping_rate/scale)*w;
        const double p=o.plasma_frequency/scale;
        const double numerator=p*p;
        if (real==0.0 && loss==0.0)
            throw std::domain_error("Optical oscillator pole or unrepresentable denominator.");
        if (!std::isfinite(numerator) || numerator==0.0
            || (o.damping_rate>0.0 && loss==0.0))
            throw std::overflow_error("Optical oscillator outside representable scaled range.");
        const auto term= numerator/std::complex<double>{real,-loss};
        epsilon+=term;
        if (!std::isfinite(epsilon.real()) || !std::isfinite(epsilon.imag()))
            throw std::overflow_error("Optical permittivity overflow.");
    }
    return epsilon;
}

// Vacuum wavelength in metres. Uses the exact SI speed of light.
[[nodiscard]] inline double optical_angular_frequency(double wavelength_m) {
    if (!(wavelength_m>0.0) || !std::isfinite(wavelength_m))
        throw std::invalid_argument("Invalid vacuum wavelength.");
    constexpr double two_pi_c=1883651567.3088532;
    const double omega=two_pi_c/wavelength_m;
    if (!(omega>0.0) || !std::isfinite(omega))
        throw std::overflow_error("Optical angular frequency outside representable range.");
    return omega;
}

} // namespace nova::physics
