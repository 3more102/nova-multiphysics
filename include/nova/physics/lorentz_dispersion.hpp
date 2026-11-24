#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <span>
#include <stdexcept>

namespace nova::optics {

// SI angular frequencies (rad/s); exp(-i omega t).
struct LorentzPole {
    double resonance_frequency;
    double damping_frequency;
    double strength; // dimensionless static susceptibility contribution
};
struct LorentzResponse {
    std::complex<double> relative_permittivity;
    std::complex<double> derivative_per_angular_frequency;
};

// Passive, isotropic, local dielectric. Lossless real-axis poles are excluded.
// The derivative is analytic, not a finite-difference approximation.
[[nodiscard]] inline LorentzResponse lorentz_permittivity(
    double omega, double epsilon_infinity, std::span<const LorentzPole> poles) {
    if (!std::isfinite(omega) || omega < 0.0 ||
        !std::isfinite(epsilon_infinity) || epsilon_infinity <= 0.0)
        throw std::invalid_argument("Invalid Lorentz frequency or background.");
    LorentzResponse out{{epsilon_infinity,0.0},{}};
    const auto finite = [](std::complex<double> z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    };
    for (const auto& p : poles) {
        if (!std::isfinite(p.resonance_frequency) || p.resonance_frequency <= 0.0 ||
            !std::isfinite(p.damping_frequency) || p.damping_frequency < 0.0 ||
            !std::isfinite(p.strength) || p.strength < 0.0)
            throw std::invalid_argument("Invalid passive Lorentz pole.");
        if (p.strength == 0.0) continue;
        if (omega == 0.0) {
            out.relative_permittivity += p.strength;
            // At DC: epsilon' = i S gamma / omega_0^2.
            const long double slope = static_cast<long double>(p.strength) *
                (static_cast<long double>(p.damping_frequency) / p.resonance_frequency) /
                p.resonance_frequency;
            out.derivative_per_angular_frequency +=
                std::complex<double>{0.0, static_cast<double>(slope)};
            continue;
        }
        if (omega == p.resonance_frequency && p.damping_frequency == 0.0)
            throw std::domain_error("Lossless Lorentz resonance is singular.");
        // Normalize before squaring to avoid overflow at optical frequencies.
        const double scale = std::max({omega,p.resonance_frequency,p.damping_frequency});
        const double w = omega/scale, r = p.resonance_frequency/scale;
        const double g = p.damping_frequency/scale;
        const std::complex<double> denominator{((p.resonance_frequency-omega)/scale)*(r+w),-g*w};
        if (denominator == std::complex<double>{})
            throw std::overflow_error("Lorentz denominator outside numerical range.");
        const auto contribution = ((r*r)/denominator)*p.strength;
        const auto derivative = (contribution/denominator)*
            std::complex<double>{2.0*w,g}/scale;
        out.relative_permittivity += contribution;
        out.derivative_per_angular_frequency += derivative;
    }
    if (!finite(out.relative_permittivity) || !finite(out.derivative_per_angular_frequency))
        throw std::overflow_error("Lorentz response outside numerical range.");
    return out;
}
} // namespace nova::optics
