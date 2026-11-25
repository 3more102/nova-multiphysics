#include "nova/physics/slab_waveguide.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace nova::physics {
namespace {
constexpr double pi = std::numbers::pi_v<double>;
bool positive(double x) { return std::isfinite(x) && x > 0.0; }
}

std::vector<SlabGuidedMode> solve_symmetric_slab(
    const SlabWaveguide& g, SlabPolarization pol) {
    if (!positive(g.vacuum_wavelength) || !positive(g.half_thickness) ||
        !positive(g.core_index) || !positive(g.cladding_index) ||
        g.core_index <= g.cladding_index ||
        (pol != SlabPolarization::te && pol != SlabPolarization::tm))
        throw std::invalid_argument("invalid symmetric dielectric slab");
    const double k0 = 2.0*pi/g.vacuum_wavelength;
    const double nc2 = g.core_index*g.core_index;
    const double ns2 = g.cladding_index*g.cladding_index;
    const double delta = (g.core_index-g.cladding_index)*
                         (g.core_index+g.cladding_index);
    const double v = (k0*g.half_thickness)*std::sqrt(delta);
    const double factor = pol == SlabPolarization::tm ? nc2/ns2 : 1.0;
    if (!positive(k0) || !positive(nc2) || !positive(ns2) ||
        !positive(delta) || !positive(factor) ||
        !std::isfinite(v) || v < 1e-4 || v > 1e4)
        throw std::invalid_argument("unsupported slab numerical scale (V in [1e-4, 1e4])");
    std::vector<SlabGuidedMode> modes;
    for (std::size_t m=0; ; ++m) {
        const double cutoff = static_cast<double>(m)*pi/2.0;
        if (cutoff >= v) break;
        // u=V cos(theta), w=V sin(theta): avoids cancellation at cutoff.
        const auto residual = [&](double theta) {
            const double s=std::sin(theta), c=std::cos(theta);
            return v*c-cutoff-std::atan2(factor*s,c);
        };
        double lo=0.0, hi=pi/2.0;
        for (int i=0; i<100; ++i) {
            const double mid=lo+(hi-lo)/2.0;
            if (mid==lo || mid==hi) break;
            if (residual(mid)>0.0) lo=mid; else hi=mid;
        }
        const double theta=lo+(hi-lo)/2.0;
        const double s=std::sin(theta), c=std::cos(theta);
        const double u=v*c, w=v*s;
        const double neff=std::sqrt(ns2+delta*s*s);
        const double beta=k0*neff, decay=w/g.half_thickness;
        const double err=std::abs(residual(theta));
        if (!positive(u) || !positive(w) || !positive(beta) ||
            !positive(decay) || !(neff>g.cladding_index && neff<g.core_index) ||
            err>1e-10*std::max(1.0,v))
            throw std::runtime_error("slab mode not representable at double precision");
        const bool even=m%2==0;
        const double sinc=std::sin(2.0*u)/(2.0*u);
        const double core=1.0+(even ? sinc : -sinc);
        const double edge=even ? std::cos(u) : std::sin(u);
        // Integrals of |E_y|^2 (TE) or |H_y|^2/epsilon (TM).
        // Common half-thickness and Poynting prefactors cancel.
        const double tail=edge*edge/w;
        const double fraction=core/(core+factor*tail);
        if (!std::isfinite(fraction) || fraction<0.0 || fraction>1.0)
            throw std::runtime_error("invalid slab power fraction");
        modes.push_back({m,pol,neff,beta,decay,u,w,fraction,err});
    }
    return modes;
}

} // namespace nova::physics
