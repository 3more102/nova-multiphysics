#include "nova/physics/shg_reference.hpp"
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {
void check(bool ok, const char* what) {
    if (!ok) throw std::runtime_error(what);
}
void near(double x, double y, double tolerance, const char* what) {
    check(std::isfinite(x) && std::abs(x-y) <= tolerance, what);
}
}
int main() {
    using namespace nova::optics;
    try {
        ShgInput p{1.064e-6, 2.1, 2.2, 10e-12, 1e8, 1e-3};
        constexpr double c = 299792458.;
        constexpr double eps0 = 8.8541878128e-12;
        const double omega = 2 * std::numbers::pi * c / p.vacuum_wavelength_m;
        const double eta = 2 * omega * omega * p.effective_d_m_per_v *
            p.effective_d_m_per_v / (p.fundamental_index * p.fundamental_index *
            p.harmonic_index * eps0 * c * c * c);
        const double expected = eta * p.pump_intensity_w_per_m2 * p.length_m * p.length_m;
        const auto u = undepleted_shg(p, 0);
        near(u.conversion_fraction, expected, expected * 1e-13, "SI normalization");
        near(u.harmonic_intensity_w_per_m2, expected * p.pump_intensity_w_per_m2,
             u.harmonic_intensity_w_per_m2 * 1e-13, "absolute harmonic intensity");
        const double mismatch = std::numbers::pi / p.length_m;
        near(undepleted_shg(p, mismatch).conversion_fraction / expected,
             4 / (std::numbers::pi * std::numbers::pi), 1e-13, "sinc envelope");
        near(undepleted_shg(p, -mismatch).conversion_fraction,
             undepleted_shg(p, mismatch).conversion_fraction, expected * 1e-13, "even mismatch");
        check(undepleted_shg(p, 2 * mismatch).conversion_fraction < expected * 1e-28,
              "first coherence null");
        auto q = p;
        q.length_m *= 2;
        near(undepleted_shg(q,0).conversion_fraction / expected, 4, 1e-13, "length scaling");
        q = p; q.pump_intensity_w_per_m2 *= 2;
        near(undepleted_shg(q,0).harmonic_intensity_w_per_m2 / u.harmonic_intensity_w_per_m2,
             4, 1e-13, "quadratic intensity scaling");
        q = p; q.effective_d_m_per_v *= -1;
        near(undepleted_shg(q,0).conversion_fraction, expected, expected * 1e-13, "sign symmetry");
        q = p; q.length_m *= 1e-3;
        const auto weak = undepleted_shg(q,0);
        const auto weak_exact = phase_matched_depleted_shg(q);
        near(weak_exact.conversion_fraction / weak.conversion_fraction, 1, 1e-7, "weak limit");
        q = p; q.length_m /= std::sqrt(expected);
        const auto exact = phase_matched_depleted_shg(q);
        near(exact.conversion_fraction, std::pow(std::tanh(1.),2), 1e-13, "depleted solution");
        near((exact.pump_intensity_w_per_m2 + exact.harmonic_intensity_w_per_m2) /
             q.pump_intensity_w_per_m2, 1, 1e-13, "energy conservation");
        q.length_m *= 20;
        const auto saturated = phase_matched_depleted_shg(q);
        check(saturated.conversion_fraction > .999999 && saturated.conversion_fraction <= 1 &&
              saturated.pump_intensity_w_per_m2 > 0, "saturation and residual pump");
        check(undepleted_shg(q,0).conversion_fraction > 1, "unclipped approximation");
        q.length_m = 0;
        near(phase_matched_depleted_shg(q).pump_intensity_w_per_m2,
             q.pump_intensity_w_per_m2, 0, "zero length");
        near(undepleted_shg(q,0).conversion_fraction, 0, 0, "zero conversion");
        q = p; q.effective_d_m_per_v = 0;
        near(phase_matched_depleted_shg(q).conversion_fraction, 0, 0, "zero nonlinearity");
        q = p; q.harmonic_index = -1;
        bool rejected = false;
        try { (void)phase_matched_depleted_shg(q); }
        catch (const std::invalid_argument&) { rejected = true; }
        check(rejected, "invalid index");
        rejected = false;
        try { (void)undepleted_shg(p, std::numeric_limits<double>::infinity()); }
        catch (const std::invalid_argument&) { rejected = true; }
        check(rejected, "invalid mismatch");
        std::cout << "NOVA SHG reference checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
