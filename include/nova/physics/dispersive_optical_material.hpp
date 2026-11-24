#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <stdexcept>

namespace nova::physics {

// All frequency parameters are angular frequencies in rad/s.
// resonance == 0 represents a Drude pole; positive resonance a Lorentz pole.
// strength is the square root of the numerator, in rad/s.
struct OpticalPole {
    double resonance{};
    double strength{};
    double damping{};
};

inline double optical_angular_frequency(double vacuum_wavelength_m) {
    if (!(vacuum_wavelength_m > 0.0) || !std::isfinite(vacuum_wavelength_m))
        throw std::invalid_argument("Invalid vacuum wavelength.");
    const double omega = 2.0 * std::numbers::pi * 299792458.0 / vacuum_wavelength_m;
    if (!(omega > 0.0) || !std::isfinite(omega))
        throw std::overflow_error("Optical angular frequency out of range.");
    return omega;
}

// exp(-i omega t): passive materials have Im(epsilon_r) >= 0.
// epsilon_r = epsilon_inf + sum strength^2/(resonance^2-omega^2-i*damping*omega).
// Scale the pole evaluation to avoid squaring dimensional optical frequencies.
inline std::complex<double> dispersive_permittivity(
    double omega, double epsilon_inf, std::span<const OpticalPole> poles) {
    if (!(omega > 0.0) || !std::isfinite(omega) ||
        !(epsilon_inf > 0.0) || !std::isfinite(epsilon_inf))
        throw std::invalid_argument("Invalid dispersive optical material or frequency.");
    std::complex<double> result{epsilon_inf, 0.0};
    for (const auto& p : poles) {
        if (!std::isfinite(p.resonance) || p.resonance < 0.0 ||
            !std::isfinite(p.strength) || p.strength < 0.0 ||
            !std::isfinite(p.damping) || p.damping < 0.0)
            throw std::invalid_argument("Invalid passive optical pole.");
        if (p.strength == 0.0) continue;
        const double scale = std::max({omega, p.resonance, p.strength, p.damping});
        const double w = omega / scale, r = p.resonance / scale;
        const double s = p.strength / scale, g = p.damping / scale;
        // Difference of squares avoids subtracting two nearly equal squares.
        const double difference = (p.resonance - omega) / scale;
        const std::complex<double> denominator{difference * (r + w), -g * w};
        if (denominator == std::complex<double>{})
            throw std::domain_error("Singular or unrepresentable optical pole.");
        if (s * s == 0.0)
            throw std::overflow_error("Optical pole strength underflow.");
        const auto contribution = (s * s) / denominator;
        result += contribution;
        if (!std::isfinite(result.real()) || !std::isfinite(result.imag()) ||
            result.imag() < 0.0)
            throw std::overflow_error("Unrepresentable passive optical permittivity.");
    }
    return result;
}

// Principal passive refractive index, for nonmagnetic isotropic materials.
// n + i*kappa gives spatial decay for a forward exp(i*k*z) wave.
inline std::complex<double> passive_refractive_index(std::complex<double> epsilon_r) {
    if (!std::isfinite(epsilon_r.real()) || !std::isfinite(epsilon_r.imag()) ||
        epsilon_r.imag() < 0.0)
        throw std::invalid_argument("Invalid passive optical permittivity.");
    // Normalize signed zero so negative real epsilon takes the +i branch.
    if (epsilon_r.imag() == 0.0) epsilon_r.imag(0.0);
    return std::sqrt(epsilon_r);
}

} // namespace nova::physics
