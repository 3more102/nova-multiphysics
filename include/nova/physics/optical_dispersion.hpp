#pragma once

#include "nova/physics/maxwell_3d.hpp"

#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <stdexcept>

namespace nova::physics {

// SI units: angular frequencies in rad/s; vacuum wavelengths in metres.
// Passive exp(-i omega t) convention, consistent with Maxwell3D.
inline constexpr double optical_speed_of_light = 299792458.0;

struct DrudePole {
    double plasma_frequency{};
    double damping{};
};

struct LorentzPole {
    double resonance_frequency{};
    double oscillator_strength{}; // dimensionless static susceptibility
    double damping{};
};

namespace optical_detail {
inline void positive(double x) {
    if (!(std::isfinite(x) && x > 0.0))
        throw std::invalid_argument("Expected finite positive optical parameter");
}
inline void nonnegative(double x) {
    if (!(std::isfinite(x) && x >= 0.0))
        throw std::invalid_argument("Expected finite nonnegative optical parameter");
}
inline std::complex<double> passive(std::complex<double> x) {
    if (!std::isfinite(x.real()) || !std::isfinite(x.imag()) || x.imag() < 0.0)
        throw std::overflow_error("Nonfinite or nonpassive optical permittivity");
    return x;
}
} // namespace optical_detail

[[nodiscard]] inline double optical_angular_frequency(double vacuum_wavelength) {
    optical_detail::positive(vacuum_wavelength);
    const double omega = (2.0 * std::numbers::pi * optical_speed_of_light) / vacuum_wavelength;
    optical_detail::positive(omega);
    return omega;
}

// epsilon = epsilon_inf - sum(wp^2/(w^2+i gamma w))
//          + sum(f w0^2/(w0^2-w^2-i gamma w)).
[[nodiscard]] inline std::complex<double> drude_lorentz_permittivity(
    double omega, double epsilon_infinity,
    std::span<const DrudePole> drude = {},
    std::span<const LorentzPole> lorentz = {}) {
    using C = std::complex<double>;
    optical_detail::positive(omega);
    optical_detail::positive(epsilon_infinity);
    C eps{epsilon_infinity, 0.0};
    for (const auto& p : drude) {
        optical_detail::nonnegative(p.plasma_frequency);
        optical_detail::nonnegative(p.damping);
        if (p.plasma_frequency == 0.0) continue;
        const double ratio = p.plasma_frequency / omega;
        const C denom{1.0, p.damping / omega};
        if (!std::isfinite(denom.imag()))
            throw std::overflow_error("Drude damping ratio outside numerical range");
        eps += optical_detail::passive(-ratio * ratio / denom);
    }
    for (const auto& p : lorentz) {
        optical_detail::positive(p.resonance_frequency);
        optical_detail::nonnegative(p.oscillator_strength);
        optical_detail::nonnegative(p.damping);
        if (p.oscillator_strength == 0.0) continue;
        const double ratio = omega / p.resonance_frequency;
        const C denom{(1.0-ratio)*(1.0+ratio), -(p.damping/p.resonance_frequency)*ratio};
        if (!std::isfinite(denom.real()) || !std::isfinite(denom.imag()))
            throw std::overflow_error("Lorentz frequency ratio outside numerical range");
        if (denom == C{})
            throw std::domain_error("Undamped Lorentz resonance is singular");
        eps += optical_detail::passive(p.oscillator_strength / denom);
    }
    return optical_detail::passive(eps);
}

// For a nonmagnetic optical medium, epsilon_r=(n+i k)^2.
// n and k are nonnegative refractive index and extinction coefficient.
[[nodiscard]] inline MaxwellCellMaterial3D optical_material_from_nk(double n, double k) {
    optical_detail::nonnegative(n);
    optical_detail::nonnegative(k);
    return {optical_detail::passive(std::complex<double>{n,k} * std::complex<double>{n,k}),1.0};
}

[[nodiscard]] inline MaxwellCellMaterial3D optical_drude_lorentz_material(
    double vacuum_wavelength, double epsilon_infinity,
    std::span<const DrudePole> drude = {},
    std::span<const LorentzPole> lorentz = {}) {
    return {drude_lorentz_permittivity(optical_angular_frequency(vacuum_wavelength),
                                       epsilon_infinity, drude, lorentz), 1.0};
}

struct OpticalIndexSample {
    double vacuum_wavelength{};
    double n{};
    double k{};
};

// Piecewise linear interpolation of n and k in vacuum wavelength.
// Require strictly ascending, finite data; never silently extrapolate.
[[nodiscard]] inline MaxwellCellMaterial3D interpolate_optical_index(
    std::span<const OpticalIndexSample> samples, double vacuum_wavelength) {
    optical_detail::positive(vacuum_wavelength);
    if (samples.empty()) throw std::invalid_argument("Empty optical index table");
    double previous = 0.0;
    for (const auto& sample : samples) {
        optical_detail::positive(sample.vacuum_wavelength);
        optical_detail::nonnegative(sample.n);
        optical_detail::nonnegative(sample.k);
        if (sample.vacuum_wavelength <= previous)
            throw std::invalid_argument("Optical wavelength table must be strictly increasing");
        previous = sample.vacuum_wavelength;
    }
    if (vacuum_wavelength < samples.front().vacuum_wavelength ||
        vacuum_wavelength > samples.back().vacuum_wavelength)
        throw std::out_of_range("Optical wavelength outside tabulated interval");
    if (vacuum_wavelength == samples.front().vacuum_wavelength)
        return optical_material_from_nk(samples.front().n, samples.front().k);
    for (std::size_t i=1; i<samples.size(); ++i) {
        const auto& a = samples[i-1];
        const auto& b = samples[i];
        if (vacuum_wavelength <= b.vacuum_wavelength) {
            const double t = (vacuum_wavelength-a.vacuum_wavelength) /
                             (b.vacuum_wavelength-a.vacuum_wavelength);
            return optical_material_from_nk(std::lerp(a.n,b.n,t),std::lerp(a.k,b.k,t));
        }
    }
    throw std::logic_error("Unreachable optical table interval");
}

} // namespace nova::physics
