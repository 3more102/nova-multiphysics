#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <stdexcept>

namespace nova::optics {
using JonesField = std::array<std::complex<double>, 2>;
using JonesMatrix = std::array<JonesField, 2>;
struct JonesLayer {
    std::complex<double> fast_index, slow_index;
    double thickness_m, axis_rad;
};
// Forward-only normal-incidence propagation; no interface reflections.
[[nodiscard]] inline JonesMatrix jones_stack(double wavelength_m, std::span<const JonesLayer> layers) {
    using C = std::complex<double>;
    if (!std::isfinite(wavelength_m) || wavelength_m <= 0)
        throw std::invalid_argument("invalid vacuum wavelength");
    const double k = 2 * std::numbers::pi_v<double> / wavelength_m;
    if (!std::isfinite(k)) throw std::overflow_error("nonfinite wavenumber");
    JonesMatrix total{{{C{1}, C{0}}, {C{0}, C{1}}}};
    for (const auto& l : layers) {
        if (!std::isfinite(l.thickness_m) || l.thickness_m < 0 || !std::isfinite(l.axis_rad))
            throw std::invalid_argument("invalid Jones layer geometry");
        const auto phase = [&](C n) {
            if (!std::isfinite(n.real()) || !std::isfinite(n.imag()) || n.real() <= 0 || n.imag() < 0)
                throw std::invalid_argument("invalid passive principal index");
            const C z = C{0,1} * n * (k * l.thickness_m);
            if (!std::isfinite(z.real()) || !std::isfinite(z.imag()))
                throw std::overflow_error("nonfinite Jones phase");
            return std::exp(z);
        };
        const C a = phase(l.fast_index), b = phase(l.slow_index);
        const double c = std::cos(l.axis_rad), s = std::sin(l.axis_rad);
        const C xy = (a-b)*c*s;
        const JonesMatrix m{{{a*c*c+b*s*s, xy}, {xy, a*s*s+b*c*c}}};
        JonesMatrix next{};
        for (std::size_t i=0;i<2;++i)
            for (std::size_t j=0;j<2;++j)
                next[i][j] = m[i][0]*total[0][j] + m[i][1]*total[1][j];
        total = next;
    }
    return total;
}
[[nodiscard]] inline JonesField apply_jones(const JonesMatrix& m, const JonesField& e) {
    return {m[0][0]*e[0]+m[0][1]*e[1], m[1][0]*e[0]+m[1][1]*e[1]};
}
} // namespace nova::optics
