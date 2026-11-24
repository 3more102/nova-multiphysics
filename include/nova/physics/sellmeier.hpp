#pragma once

#include <cmath>
#include <numbers>
#include <span>
#include <stdexcept>

namespace nova::optics {

// Sellmeier coefficients B are dimensionless; resonance wavelengths are in metres.
struct SellmeierTerm {
    double strength;
    double resonance_wavelength_m;
};

struct SellmeierResult {
    double phase_index;
    double group_index;
    double dn_dlambda_per_m;
    double d2n_dlambda2_per_m2;
    double beta2_s2_per_m;
    double dispersion_s_per_m2;
};

// Transparent isotropic material reference model (real n > 0).
// The caller is responsible for staying within the fit's transparent validity band.
// No absorption, anisotropy, gain, or waveguide/geometrical dispersion is modeled.
[[nodiscard]] inline SellmeierResult sellmeier_dispersion(
    double wavelength_m, std::span<const SellmeierTerm> terms) {
    if (!std::isfinite(wavelength_m) || wavelength_m <= 0.0)
        throw std::invalid_argument("invalid vacuum wavelength");
    double f = 1.0, f1 = 0.0, f2 = 0.0;
    // Accumulate f=n^2 and wavelength-scaled derivatives to avoid
    // unnecessarily large intermediate inverse powers of wavelength.
    for (const auto& term : terms) {
        if (!std::isfinite(term.strength) ||
            !std::isfinite(term.resonance_wavelength_m) ||
            term.resonance_wavelength_m <= 0.0)
            throw std::invalid_argument("invalid Sellmeier term");
        if (term.strength == 0.0) continue;
        const double ratio = term.resonance_wavelength_m / wavelength_m;
        const double u = ratio * ratio;
        if (!std::isfinite(u))
            throw std::overflow_error("Sellmeier resonance ratio overflow");
        const double d = 1.0 - u;
        if (d == 0.0) throw std::domain_error("Sellmeier resonance pole");
        const double v = term.strength / d;
        const double w = u / d;
        f += v;
        f1 -= 2.0 * v * w;
        f2 += 2.0 * v * w * ((3.0 + u) / d);
    }
    if (!std::isfinite(f) || !std::isfinite(f1) || !std::isfinite(f2))
        throw std::overflow_error("nonfinite Sellmeier response");
    if (f <= 0.0) throw std::domain_error("nonpositive real squared index");
    const double n = std::sqrt(f);
    const double n1 = 0.5 * f1 / n;
    const double n2 = 0.5 * f2 / n - n1 * (n1 / n);
    constexpr double c = 299792458.0;
    const SellmeierResult result{
        n, n - n1, n1 / wavelength_m,
        (n2 / wavelength_m) / wavelength_m,
        (wavelength_m / (2.0 * std::numbers::pi_v<double> * c * c)) * n2,
        -(n2 / c) / wavelength_m
    };
    if (!std::isfinite(result.group_index) ||
        !std::isfinite(result.dn_dlambda_per_m) ||
        !std::isfinite(result.d2n_dlambda2_per_m2) ||
        !std::isfinite(result.beta2_s2_per_m) ||
        !std::isfinite(result.dispersion_s_per_m2))
        throw std::overflow_error("nonfinite Sellmeier dispersion");
    return result;
}

} // namespace nova::optics
