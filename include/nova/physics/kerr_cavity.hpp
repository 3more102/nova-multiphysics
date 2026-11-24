#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace nova::optics {
// Single-mode normalized dynamics: da/dtau = -(1+i*(detuning-|a|^2))*a + sqrt(drive).
// x=|a|^2 satisfies drive=x*(1+(detuning-x)^2).
enum class KerrCavityStability { stable, unstable, marginal };
struct KerrCavityState {
    double intensity;
    double differential_slope;
    KerrCavityStability stability;
};
// All distinct nonnegative steady states, in ascending intensity order.
// Bounded reference model, not a multimode or spatial nonlinear solver.
[[nodiscard]] inline std::vector<KerrCavityState> kerr_cavity_steady_states(
    double detuning, double drive) {
    if (!std::isfinite(detuning) || !std::isfinite(drive) ||
        std::abs(detuning) > 1000.0 || drive < 0.0 || drive > 1e6 ||
        (drive > 0.0 && drive < 1e-12))
        throw std::invalid_argument("Kerr cavity input outside reference range");
    const long double d = detuning, p = drive;
    const auto residual = [=](long double x) {
        const auto z = d-x;
        return x*(1+z*z)-p;
    };
    std::vector<KerrCavityState> states;
    const auto append = [&](long double x) {
        const long double slope = 1+(d-x)*(d-3*x);
        const long double scale = 1+std::abs((d-x)*(d-3*x));
        const long double tol = 64*std::numeric_limits<double>::epsilon()*scale;
        const auto stability = std::abs(slope) <= tol ? KerrCavityStability::marginal :
            (slope > 0 ? KerrCavityStability::stable : KerrCavityStability::unstable);
        const double value = static_cast<double>(x);
        if (states.empty() || states.back().intensity != value)
            states.push_back({value, static_cast<double>(slope), stability});
    };
    if (p == 0) { append(0); return states; }
    // x<=p; split into monotone intervals at the two turning points.
    std::vector<long double> cuts{0};
    if (d > 0 && d*d > 3) {
        const long double q = std::sqrt(d*d-3);
        for (const auto x : {(2*d-q)/3, (2*d+q)/3})
            if (x > 0 && x < p) cuts.push_back(x);
    }
    cuts.push_back(p);
    for (std::size_t i=1; i<cuts.size(); ++i) {
        auto lo=cuts[i-1], hi=cuts[i];
        auto fl=residual(lo), fh=residual(hi);
        if (fl == 0) append(lo);
        if ((fl < 0 && fh > 0) || (fl > 0 && fh < 0)) {
            for (int step=0; step<160; ++step) {
                const auto mid=std::midpoint(lo,hi);
                if (mid == lo || mid == hi) break;
                const auto fm=residual(mid);
                if (fm == 0) { lo=hi=mid; break; }
                if ((fl < 0) == (fm < 0)) { lo=mid; fl=fm; }
                else hi=mid;
            }
            append(std::midpoint(lo,hi));
        }
        if (fh == 0) append(cuts[i]);
    }
    return states;
}
} // namespace nova::optics
