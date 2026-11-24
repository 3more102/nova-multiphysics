#pragma once
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace nova::optics {
enum class SlabPolarization { TE, TM };

struct SlabMode {
    std::size_t order;
    double effective_index;
    double beta_per_m;
    double core_transverse_wavenumber_per_m;
    double cladding_decay_per_m;
    double normalized_dispersion_residual;
};

// Symmetric, lossless, nonmagnetic step-index slab, infinite claddings.
// thickness_m is the FULL core thickness. Vacuum wavelength and SI units.
// Returns bound modes in descending effective-index order; cutoff excluded.
[[nodiscard]] inline std::vector<SlabMode> symmetric_slab_modes(
    double wavelength_m, double thickness_m, double core_index,
    double cladding_index, SlabPolarization polarization) {
    constexpr double pi = std::numbers::pi_v<double>;
    auto positive = [](double x) { return std::isfinite(x) && x > 0.0; };
    if (!positive(wavelength_m) || !positive(thickness_m) ||
        !positive(core_index) || !positive(cladding_index) ||
        core_index <= cladding_index ||
        (polarization != SlabPolarization::TE && polarization != SlabPolarization::TM))
        throw std::invalid_argument("invalid symmetric slab parameters");
    const double a = thickness_m / 2.0;
    const double k0 = 2.0 * pi / wavelength_m;
    const double contrast = (core_index-cladding_index)*(core_index+cladding_index);
    const double v = (k0*a)*std::sqrt(contrast);
    const double ratio = core_index/cladding_index;
    const double q = polarization == SlabPolarization::TE ? 1.0 : ratio*ratio;
    if (!positive(a) || !positive(k0) || !positive(contrast) ||
        !positive(v) || !positive(q) || v > 10000.0)
        throw std::domain_error("slab outside supported numerical range");
    std::vector<SlabMode> result;
    for (std::size_t m=0; ; ++m) {
        const double cutoff = static_cast<double>(m)*pi/2.0;
        if (!(v > cutoff)) break;
        // u=V cos(theta), w=V sin(theta), so no cancellation near cutoff.
        const auto f = [&](double theta) {
            return v*std::cos(theta)-cutoff-
                std::atan2(q*std::sin(theta), std::cos(theta));
        };
        double lo=0.0, hi=pi/2.0;
        for (int i=0; i<100; ++i) {
            const double mid=lo+(hi-lo)/2.0;
            if (mid==lo || mid==hi) break;
            if (f(mid)>0.0) lo=mid; else hi=mid;
        }
        const double theta=lo+(hi-lo)/2.0;
        const double h=(v*std::cos(theta))/a;
        const double decay=(v*std::sin(theta))/a;
        const double neff=std::hypot(cladding_index, decay/k0);
        const double beta=k0*neff;
        const double residual=std::abs(f(theta));
        if (!positive(h) || !positive(decay) || !positive(beta) ||
            !(neff>cladding_index && neff<core_index) ||
            !std::isfinite(residual) || residual>1e-10)
            throw std::domain_error("slab mode unresolved at floating-point precision");
        result.push_back({m, neff, beta, h, decay, residual});
    }
    return result;
}
} // namespace nova::optics
