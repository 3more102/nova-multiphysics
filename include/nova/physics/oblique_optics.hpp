#pragma once
#include "nova/physics/optics.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace nova::optics {
enum class Polarization { TE, TM };
// Tangential electric-field amplitudes; angle in radians, SI units.
// Exact zero longitudinal wavenumber is outside the supported domain.
[[nodiscard]] inline StackResult oblique_incidence_stack(
    double wavelength, double incident, double exit, double angle,
    Polarization polarization, std::span<const Layer> layers) {
    const auto finite = [](Complex z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    };
    const auto checked = [&](Complex z) {
        if (!finite(z)) throw std::overflow_error("nonfinite optical result");
        return z;
    };
    const auto positive = [](double x) { return std::isfinite(x) && x > 0; };
    if (!positive(wavelength) || !positive(incident) || !positive(exit) ||
        !std::isfinite(angle) || angle < 0 || angle >= std::numbers::pi_v<double>/2 ||
        (polarization != Polarization::TE && polarization != Polarization::TM))
        throw std::invalid_argument("invalid oblique incidence parameters");
    const double k = 2*std::numbers::pi_v<double>/wavelength;
    const double s = incident*std::sin(angle);
    if (!positive(k) || !std::isfinite(s*s)) throw std::overflow_error("invalid optical scale");
    const auto longitudinal = [&](Complex n) {
        const Complex v = checked(n*n-s*s);
        const Complex q = std::sqrt(Complex{v.real(), std::abs(v.imag())});
        if (q == Complex{}) throw std::domain_error("zero longitudinal wavenumber");
        return checked(q);
    };
    const auto admittance = [&](Complex n, Complex q) {
        return checked(polarization == Polarization::TE ? q : n*n/q);
    };
    const auto interface = [&](Complex a, Complex b) {
        const Complex sum = checked(a+b);
        if (sum == Complex{}) throw std::domain_error("singular interface");
        return checked((a-b)/sum);
    };
    const Complex yi = admittance(incident, longitudinal(incident));
    const Complex ye = admittance(exit, longitudinal(exit));
    if (!positive(yi.real())) throw std::domain_error("invalid incident flux");
    Complex right = ye, r = 0.0, t = 1.0;
    for (std::size_t j = layers.size(); j > 0; --j) {
        const auto& l = layers[j-1];
        if (!finite(l.index) || l.index.real() < 0 || l.index.imag() < 0 ||
            l.index == Complex{} || !std::isfinite(l.thickness_m) || l.thickness_m < 0)
            throw std::invalid_argument("invalid passive layer");
        if (l.thickness_m == 0) continue;
        const Complex q = longitudinal(l.index);
        const Complex y = admittance(l.index,q);
        const Complex ri = interface(y,right);
        const Complex den = checked(1.0+ri*r);
        if (den == Complex{}) throw std::domain_error("singular stack");
        const Complex p = std::exp(checked(Complex{0,1}*q*(k*l.thickness_m)));
        t = checked(p*(1.0+ri)*t/den);
        r = checked(p*p*(ri+r)/den);
        right = y;
    }
    const Complex ri = interface(yi,right);
    const Complex den = checked(1.0+ri*r);
    if (den == Complex{}) throw std::domain_error("singular entrance");
    t = checked((1.0+ri)*t/den);
    r = checked((ri+r)/den);
    const double R = std::norm(r), T = ye.real()/yi.real()*std::norm(t);
    const double A = 1-R-T;
    if (!std::isfinite(R) || !std::isfinite(T) || !std::isfinite(A))
        throw std::overflow_error("nonfinite optical power");
    if (A < -1e-9 || T < 0) throw std::runtime_error("passive optical power balance failed");
    return {r,t,R,T,A};
}
} // namespace nova::optics
