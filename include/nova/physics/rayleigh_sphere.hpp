#pragma once
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace nova::optics {

// Efficiencies are cross sections divided by pi * radius^2.
// The dielectric contrast is epsilon_particle / epsilon_host.
struct RayleighSphereResult {
    double size_parameter;
    std::complex<double> dipole_factor;
    double scattering_efficiency;
    double absorption_efficiency;
    double extinction_efficiency;
    double scattering_cross_section_m2;
    double absorption_cross_section_m2;
    double extinction_cross_section_m2;
};

// Radiatively corrected electrostatic electric dipole; exp(-i omega t).
// Lossless, nonmagnetic host. This is not a full Mie solution.
// Restricted to x <= 0.1 and |m|*x <= 0.1; near-pole cases are rejected.
[[nodiscard]] inline RayleighSphereResult rayleigh_sphere(
    double vacuum_wavelength_m, double radius_m, double host_index,
    std::complex<double> permittivity_contrast) {
    using C = std::complex<double>;
    const auto positive = [](double v) { return std::isfinite(v) && v > 0.0; };
    const auto finite = [](C v) {
        return std::isfinite(v.real()) && std::isfinite(v.imag());
    };
    if (!positive(vacuum_wavelength_m) || !positive(radius_m) ||
        !positive(host_index) || !finite(permittivity_contrast) ||
        permittivity_contrast.imag() < 0.0)
        throw std::invalid_argument("invalid passive sphere parameters");

    constexpr double pi = std::numbers::pi_v<double>;
    const double x = (2.0 * pi) * (radius_m / vacuum_wavelength_m) * host_index;
    const double internal_x = x * std::abs(std::sqrt(permittivity_contrast));
    if (!positive(x) || x > 0.1 || !std::isfinite(internal_x) || internal_x > 0.1)
        throw std::domain_error("outside restricted electric-dipole regime");

    const C numerator = permittivity_contrast - 1.0;
    const C static_den = permittivity_contrast + 2.0;
    const double scale = std::fmax(1.0, std::abs(permittivity_contrast));
    // Conservative guard against omitted O(x^2) dynamic depolarization.
    // This is a scope restriction, not an error bound.
    if (std::abs(static_den) <= 10.0 * x * x * scale)
        throw std::domain_error("near electrostatic pole; full Mie model required");

    const double g = (2.0 / 3.0) * x * x * x;
    const C den = static_den - C{0.0, g} * numerator;
    const C f = numerator / den;
    const double qsca = (8.0 / 3.0) * x * x * x * x * std::norm(f);
    // Evaluate material absorption directly, avoiding subtraction cancellation.
    const double qabs = 12.0 * x * permittivity_contrast.imag() / std::norm(den);
    const double qext = qsca + qabs;
    const double area = pi * radius_m * radius_m;
    const double sca = area * qsca, abs = area * qabs, ext = area * qext;
    if (!positive(area) || !finite(f) || !std::isfinite(qext) ||
        !std::isfinite(sca) || !std::isfinite(abs) || !std::isfinite(ext))
        throw std::overflow_error("unrepresentable sphere response");
    return {x, f, qsca, qabs, qext, sca, abs, ext};
}
} // namespace nova::optics
