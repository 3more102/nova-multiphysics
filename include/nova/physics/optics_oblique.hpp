#pragma once
#include "nova/physics/optics.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace nova::optics {
enum class Polarization { te, tm };

// Amplitudes are tangential E for TE and tangential H for TM.
// Lossless positive exterior indices; passive, nonzero layer indices.
// The exit half-space may be evanescent. Exact internal cutoff is excluded.
[[nodiscard]] inline StackResult oblique_incidence_stack(
    double wavelength_m, double incident_index, double exit_index,
    double angle_rad, Polarization polarization, std::span<const Layer> layers) {
    const auto finite = [](Complex z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    };
    if (!(std::isfinite(wavelength_m) && wavelength_m > 0 &&
          std::isfinite(incident_index) && incident_index > 0 &&
          std::isfinite(exit_index) && exit_index > 0 &&
          std::isfinite(angle_rad) && angle_rad >= 0 &&
          angle_rad < std::numbers::pi_v<double>/2) ||
        (polarization != Polarization::te && polarization != Polarization::tm))
        throw std::invalid_argument("invalid oblique stack input");
    const double k = 2*std::numbers::pi_v<double>/wavelength_m;
    const double s = incident_index*std::sin(angle_rad);
    struct Medium { Complex y, phase; };
    const auto medium = [&](Complex n, double d) {
        const auto eps = n*n;
        if (!finite(n) || n.real()<0 || n.imag()<0 || !finite(eps) ||
            eps == Complex{} || !std::isfinite(d) || d<0)
            throw std::invalid_argument("invalid passive layer");
        const auto q = std::sqrt(eps-Complex{s*s,0});
        if (!finite(q) || q == Complex{})
            throw std::invalid_argument("unsupported exact longitudinal cutoff");
        const auto y = polarization == Polarization::te ? q : q/eps;
        const auto phase = std::exp(Complex{0,1}*k*q*d);
        if (!finite(y) || !finite(phase) || y == Complex{})
            throw std::invalid_argument("unresolved optical admittance or phase");
        return Medium{y,phase};
    };
    std::vector<Medium> m;
    m.reserve(layers.size()+2);
    m.push_back(medium(incident_index,0));
    for (const auto& layer: layers) {
        if (!finite(layer.index) || layer.index.real()<0 || layer.index.imag()<0 ||
            layer.index == Complex{} || !std::isfinite(layer.thickness_m) ||
            layer.thickness_m<0)
            throw std::invalid_argument("invalid passive layer");
        if (layer.thickness_m>0) m.push_back(medium(layer.index,layer.thickness_m));
    }
    m.push_back(medium(exit_index,0));
    Complex r{}, t{1,0};
    for (std::size_t j=m.size()-1;j>0;--j) {
        const auto& a=m[j-1]; const auto& b=m[j];
        const auto interface_den=a.y+b.y;
        if (interface_den==Complex{}) throw std::runtime_error("singular optical interface");
        const auto ri=(a.y-b.y)/interface_den;
        const auto ti=2.0*a.y/interface_den;
        const auto downstream=r*b.phase*b.phase;
        const auto den=1.0+ri*downstream;
        if (den==Complex{}) throw std::runtime_error("singular optical cascade");
        t=ti*b.phase*t/den;
        r=(ri+downstream)/den;
    }
    const double R=std::norm(r);
    const double T=(m.back().y.real()/m.front().y.real())*std::norm(t);
    if (!finite(r) || !finite(t) || !std::isfinite(R) || !std::isfinite(T))
        throw std::runtime_error("nonfinite optical scattering");
    return {r,t,R,T,1-R-T};
}
} // namespace nova::optics
