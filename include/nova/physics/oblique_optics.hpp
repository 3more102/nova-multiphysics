#pragma once

#include "nova/physics/optics.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace nova::optics {

enum class Polarization { TE, TM };

// Isotropic, nonmagnetic planar layers; exp(-i omega t).
// Amplitudes are tangential electric-field ratios at the stack boundaries.
// angle_rad is measured from the normal in the lossless incident medium.
// Evanescent exit waves have zero transmitted normal power.
// Exact grazing/critical longitudinal wavenumbers and singular interfaces
// are deliberately rejected rather than silently regularized.
[[nodiscard]] inline StackResult oblique_incidence_stack(
    double wavelength_m, double incident_index, double exit_index,
    double angle_rad, Polarization polarization, std::span<const Layer> layers) {
    const auto finite = [](Complex z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    };
    const auto positive = [](double x) { return std::isfinite(x) && x > 0.0; };
    if (!positive(wavelength_m) || !positive(incident_index) ||
        !positive(exit_index) || !std::isfinite(angle_rad) || angle_rad < 0.0 ||
        angle_rad >= std::numbers::pi_v<double>/2.0 ||
        (polarization != Polarization::TE && polarization != Polarization::TM))
        throw std::invalid_argument("invalid oblique optical stack input");
    const double k = 2.0 * std::numbers::pi_v<double> / wavelength_m;
    const double s = incident_index * std::sin(angle_rad);
    const double s2 = s*s;
    if (!positive(k) || !std::isfinite(s2))
        throw std::overflow_error("nonfinite optical wavevector");
    struct Wave { Complex q; Complex y; };
    const auto wave = [&](Complex n) -> Wave {
        const Complex eps = n*n;
        const Complex rad = eps - s2;
        if (!finite(rad)) throw std::overflow_error("nonfinite longitudinal dispersion");
        // Normalize signed zero: the outgoing passive square root has Im(q)>=0.
        const Complex q = std::sqrt(Complex{rad.real(), std::abs(rad.imag())});
        if (q == Complex{}) throw std::domain_error("exact critical longitudinal wavevector");
        const Complex y = polarization == Polarization::TE ? q : eps/q;
        if (!finite(y) || y == Complex{})
            throw std::domain_error("invalid optical admittance");
        return {q,y};
    };
    for (const auto& layer : layers) {
        if (!finite(layer.index) || layer.index.real() < 0.0 ||
            layer.index.imag() < 0.0 || layer.index == Complex{} ||
            !std::isfinite(layer.thickness_m) || layer.thickness_m < 0.0)
            throw std::invalid_argument("invalid passive optical layer");
    }
    const auto entrance = wave(incident_index);
    const auto exit = wave(exit_index);
    if (!(entrance.y.real() > 0.0))
        throw std::domain_error("incident wave carries no normal power");

    Complex right = exit.y, r = 0.0, t = 1.0;
    const auto interface = [&](Complex left) {
        const Complex sum = left + right;
        if (!finite(sum) || sum == Complex{})
            throw std::domain_error("singular optical interface");
        const Complex ri = (left-right)/sum;
        const Complex den = 1.0 + ri*r;
        if (!finite(den) || den == Complex{})
            throw std::domain_error("singular optical cascade");
        // The explicit 2*left/sum avoids cancellation in 1+ri.
        t = (2.0*(left/sum))*t/den;
        r = (ri+r)/den;
    };
    for (std::size_t j=layers.size(); j>0; --j) {
        const auto& layer = layers[j-1];
        if (layer.thickness_m == 0.0) continue;
        const auto w = wave(layer.index);
        interface(w.y);
        const Complex phase = Complex{0.0,1.0}*w.q*(k*layer.thickness_m);
        if (!finite(phase)) throw std::overflow_error("nonfinite propagation phase");
        const Complex p = std::exp(phase);
        r *= p*p;
        t *= p;
        right = w.y;
        if (!finite(r) || !finite(t))
            throw std::overflow_error("nonfinite optical scattering");
    }
    interface(entrance.y);
    if (!finite(r) || !finite(t))
        throw std::overflow_error("nonfinite optical scattering");
    const double R = std::norm(r);
    const double T = exit.q.real() == 0.0 ? 0.0 :
        (exit.y.real()/entrance.y.real())*std::norm(t);
    const double A = 1.0-R-T;
    if (!std::isfinite(R) || !std::isfinite(T) || !std::isfinite(A) ||
        T < 0.0 || A < -1e-9)
        throw std::runtime_error("passive oblique stack power balance failed");
    return {r,t,R,T,A};
}

} // namespace nova::optics
