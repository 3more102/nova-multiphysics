#include "nova/physics/optics.hpp"
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace nova::optics {
namespace {
constexpr double two_pi = 2.0 * std::numbers::pi_v<double>;
bool finite(Complex x) { return std::isfinite(x.real()) && std::isfinite(x.imag()); }
void positive(double x) {
    if (!std::isfinite(x) || x <= 0.0) throw std::invalid_argument("expected finite positive input");
}
void nonnegative(double x) {
    if (!std::isfinite(x) || x < 0.0) throw std::invalid_argument("expected finite nonnegative input");
}
void passive(Complex z) {
    if (!finite(z) || z.imag() < 0.0) throw std::invalid_argument("expected finite passive permittivity");
}
double wavenumber(double wavelength) {
    positive(wavelength);
    const double k = two_pi / wavelength;
    positive(k);
    return k;
}
void check(Complex z) {
    if (!finite(z)) throw std::overflow_error("nonfinite optical result");
}
Complex interface_r(Complex a, Complex b) {
    if (a + b == Complex{}) throw std::domain_error("singular interface");
    return (a - b) / (a + b);
}
}
Complex passive_index(Complex eps) {
    passive(eps);
    // Normalize signed zero for the passive branch on the negative real axis.
    return std::sqrt(Complex{eps.real(), std::abs(eps.imag())});
}
Complex drude_permittivity(double omega, double eps_inf, double wp, double gamma) {
    positive(omega); positive(eps_inf); nonnegative(wp); nonnegative(gamma);
    const double p = wp / omega;
    const Complex eps = eps_inf - (p * p) / Complex{1.0, gamma / omega};
    check(eps);
    return eps;
}
StackResult normal_incidence_stack(double wavelength, double incident, double exit,
                                  std::span<const Layer> layers) {
    const double k = wavenumber(wavelength);
    positive(incident); positive(exit);
    for (const auto& l : layers) {
        if (!finite(l.index) || l.index.real() < 0.0 || l.index.imag() < 0.0 ||
            l.index == Complex{}) throw std::invalid_argument("invalid passive layer index");
        nonnegative(l.thickness_m);
    }
    Complex right = exit;
    Complex r = 0.0, t = 1.0;
    // Cascade from the exit backwards. Only decaying propagation factors are used.
    for (std::size_t j = layers.size(); j > 0; --j) {
        const auto& l = layers[j - 1];
        // A zero-width layer has no propagation or physical interface pair.
        // Skip before computing interface ratios to avoid catastrophic cancellation.
        if (l.thickness_m == 0.0) continue;
        const Complex ri = interface_r(l.index, right);
        const Complex den = 1.0 + ri * r;
        if (den == Complex{}) throw std::domain_error("singular stack");
        const Complex phase = Complex{0.0, 1.0} * l.index * (k * l.thickness_m);
        check(phase);
        const Complex p = std::exp(phase);
        t = p * (1.0 + ri) * t / den;
        r = p * p * (ri + r) / den;
        check(r); check(t);
        right = l.index;
    }
    const Complex ri = interface_r(incident, right);
    const Complex den = 1.0 + ri * r;
    if (den == Complex{}) throw std::domain_error("singular entrance interface");
    t = (1.0 + ri) * t / den;
    r = (ri + r) / den;
    check(r); check(t);
    const double R = std::norm(r), T = (exit / incident) * std::norm(t);
    if (!std::isfinite(R) || !std::isfinite(T)) throw std::overflow_error("nonfinite optical power");
    const double A = 1.0 - R - T;
    if (A < -1e-9) throw std::runtime_error("passive stack power balance failed");
    return {r, t, R, T, A};
}
KerrResult kerr_phase(double wavelength, double n2, double intensity, double length) {
    const double k = wavenumber(wavelength);
    if (!std::isfinite(n2)) throw std::invalid_argument("nonfinite Kerr coefficient");
    nonnegative(intensity); nonnegative(length);
    const double dn = n2 * intensity, phase = k * dn * length;
    if (!std::isfinite(dn) || !std::isfinite(phase)) throw std::overflow_error("nonfinite Kerr response");
    return {dn, phase};
}
SurfacePlasmon planar_surface_plasmon(double wavelength, double ed, Complex em) {
    const double k = wavenumber(wavelength);
    positive(ed); passive(em);
    if (em.real() >= -ed) throw std::invalid_argument("outside supported bound SPP regime");
    const Complex b = std::sqrt(ed * (em / (ed + em)));
    const Complex kd = std::sqrt(b * b - ed), km = std::sqrt(b * b - em);
    const Complex beta = k * b, decay_d = k * kd, decay_m = k * km;
    check(beta); check(decay_d); check(decay_m);
    if (beta.real() <= 0.0 || beta.imag() < 0.0 ||
        decay_d.real() <= 0.0 || decay_m.real() <= 0.0)
        throw std::domain_error("unresolved bound SPP branch");
    const double length = beta.imag() == 0.0 ? std::numeric_limits<double>::infinity()
        : 0.5 / beta.imag();
    return {beta, decay_d, decay_m, length};
}

StackResult oblique_incidence_stack(double wavelength, double incident, double exit,
                                   double angle, Polarization pol,
                                   std::span<const Layer> layers) {
    const double k = wavenumber(wavelength);
    positive(incident); positive(exit); nonnegative(angle);
    if (angle >= std::numbers::pi_v<double>/2 ||
        (pol != Polarization::te && pol != Polarization::tm))
        throw std::invalid_argument("invalid incidence angle or polarization");
    const double tangential = incident * std::sin(angle);
    const double tangential2 = tangential * tangential;
    if (!std::isfinite(tangential2)) throw std::overflow_error("nonfinite tangential wavevector");
    const auto admittance = [&](Complex n) {
        const Complex eps = n*n;
        check(eps);
        const Complex radicand = eps - tangential2;
        // Positive real or positive imaginary longitudinal branch.
        const Complex q = std::sqrt(Complex{radicand.real(), std::abs(radicand.imag())});
        check(q);
        if (q == Complex{}) throw std::domain_error("exact critical angle unsupported");
        const Complex y = pol == Polarization::te ? q : eps/q;
        check(y);
        return std::pair{q, y};
    };
    const auto [qi, yi] = admittance(incident);
    const auto [qe, ye] = admittance(exit);
    (void)qi; (void)qe;
    if (yi.real() <= 0.0) throw std::domain_error("unresolved incident power");
    Complex right = ye, r = 0.0, t = 1.0;
    for (std::size_t j = layers.size(); j > 0; --j) {
        const auto& l = layers[j-1];
        if (!finite(l.index) || l.index.real() < 0.0 || l.index.imag() < 0.0 ||
            l.index == Complex{}) throw std::invalid_argument("invalid passive layer index");
        nonnegative(l.thickness_m);
        if (l.thickness_m == 0.0) continue;
        const auto [q, y] = admittance(l.index);
        const Complex ri = interface_r(y, right);
        const Complex den = 1.0 + ri*r;
        if (den == Complex{}) throw std::domain_error("singular stack");
        const Complex phase = Complex{0.0,1.0}*q*(k*l.thickness_m);
        check(phase);
        const Complex p = std::exp(phase);
        t = p*(1.0+ri)*t/den;
        r = p*p*(ri+r)/den;
        check(r); check(t);
        right = y;
    }
    const Complex ri = interface_r(yi, right);
    const Complex den = 1.0 + ri*r;
    if (den == Complex{}) throw std::domain_error("singular entrance interface");
    t = (1.0+ri)*t/den;
    r = (ri+r)/den;
    check(r); check(t);
    const double R = std::norm(r);
    const double T = (ye.real()/yi.real())*std::norm(t);
    const double A = 1.0-R-T;
    if (!std::isfinite(R) || !std::isfinite(T) || T < 0.0)
        throw std::overflow_error("invalid optical power");
    if (A < -1e-9) throw std::runtime_error("passive stack power balance failed");
    return {r,t,R,T,A};
}
} // namespace nova::optics
