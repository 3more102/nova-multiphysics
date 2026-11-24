#include "nova/physics/sellmeier.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
void near(double a, double b, double tol, const char* message) {
    require(std::isfinite(a) && std::abs(a-b) <= tol, message);
}
template<class Exception, class Callable>
void rejects(Callable fn, const char* message) {
    try { fn(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}
}

int main() {
    using namespace nova::optics;
    try {
        const auto vacuum = sellmeier_dispersion(1.55e-6, {});
        near(vacuum.phase_index, 1.0, 0.0, "vacuum n");
        near(vacuum.group_index, 1.0, 0.0, "vacuum ng");
        near(vacuum.beta2_s2_per_m, 0.0, 0.0, "vacuum GVD");
        near(vacuum.dispersion_s_per_m2, 0.0, 0.0, "vacuum D");

        // Malitson fused silica fit: resonance wavelengths in metres.
        const std::array<SellmeierTerm, 3> silica{{
            {0.6961663, 0.0684043e-6},
            {0.4079426, 0.1162414e-6},
            {0.8974794, 9.896161e-6}
        }};
        constexpr double l = 1.55e-6, h = 1e-9;
        const auto s = sellmeier_dispersion(l, silica);
        near(s.phase_index, 1.44402362, 1e-8, "silica n");
        near(s.group_index, 1.462596, 2e-6, "silica ng");
        require(s.beta2_s2_per_m < 0.0, "silica anomalous material GVD");
        require(s.dispersion_s_per_m2 * 1e6 > 21.0 &&
                s.dispersion_s_per_m2 * 1e6 < 23.0, "silica D in ps/(nm km)");
        const double lo = sellmeier_dispersion(l-h, silica).phase_index;
        const double hi = sellmeier_dispersion(l+h, silica).phase_index;
        near(s.dn_dlambda_per_m, (hi-lo)/(2*h), 0.1, "first derivative");
        near(s.d2n_dlambda2_per_m2, (hi-2*s.phase_index+lo)/(h*h),
             1e5, "second derivative");
        constexpr double c = 299792458.0;
        near(s.dispersion_s_per_m2,
             -2*std::numbers::pi_v<double>*c*s.beta2_s2_per_m/(l*l),
             1e-17, "D and beta2 identity");
        const std::array<SellmeierTerm, 1> pole{{{1.0, l}}};
        rejects<std::domain_error>([&]{ (void)sellmeier_dispersion(l,pole); }, "pole");
        const std::array<SellmeierTerm, 1> negative{{{-3.0,l/2}}};
        rejects<std::domain_error>([&]{ (void)sellmeier_dispersion(l,negative); }, "negative n squared");
        rejects<std::invalid_argument>([&]{ (void)sellmeier_dispersion(0,silica); }, "zero wavelength");
        rejects<std::invalid_argument>([&]{ (void)sellmeier_dispersion(
            std::numeric_limits<double>::quiet_NaN(),silica); }, "nan wavelength");
        const std::array<SellmeierTerm, 1> bad{{{1.0,-l}}};
        rejects<std::invalid_argument>([&]{ (void)sellmeier_dispersion(l,bad); }, "negative resonance");
        std::cout << "NOVA Sellmeier dispersion reference cases passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
