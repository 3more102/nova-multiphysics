#include "nova/physics/rayleigh_sphere.hpp"
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {
void require(bool b, const char* message) {
    if (!b) throw std::runtime_error(message);
}
void near(double actual, double expected, double tol, const char* message) {
    require(std::isfinite(actual) && std::abs(actual - expected) <= tol, message);
}
template<class F> void rejected(F f, const char* message) {
    bool caught = false;
    try { f(); } catch (const std::invalid_argument&) { caught = true; }
    catch (const std::domain_error&) { caught = true; }
    require(caught, message);
}
}

int main() {
    using nova::optics::rayleigh_sphere;
    constexpr double lambda = 1e-6, radius = 5e-9;
    constexpr double pi = std::numbers::pi_v<double>;
    try {
        const auto matched = rayleigh_sphere(lambda, radius, 1.5, {1,0});
        near(matched.extinction_efficiency, 0, 0, "index-matched sphere");
        near(matched.size_parameter, 2*pi*1.5*radius/lambda, 1e-15, "host wavenumber");

        const auto d = rayleigh_sphere(lambda, radius, 1, {2.25,0});
        const double x = d.size_parameter, f0 = 1.25/4.25;
        const double q0 = (8.0/3.0)*std::pow(x,4)*f0*f0;
        const double g = (2.0/3.0)*std::pow(x,3);
        near(d.scattering_efficiency, q0/(1+g*g*f0*f0), 1e-16, "dipole scattering");
        near(d.absorption_efficiency, 0, 0, "lossless absorption");
        near(d.extinction_efficiency, 4*x*d.dipole_factor.imag(), 1e-16, "lossless optical theorem");

        const auto twice = rayleigh_sphere(lambda, 2*radius, 1, {2.25,0});
        near(twice.scattering_cross_section_m2/d.scattering_cross_section_m2, 64, 1e-5, "radius^6 limit");

        const auto metal = rayleigh_sphere(lambda, radius, 1, {-4,0.3});
        require(metal.absorption_efficiency > 0 && metal.scattering_efficiency > 0, "passive lossy metal");
        near(metal.extinction_efficiency, 4*metal.size_parameter*metal.dipole_factor.imag(), 1e-14, "lossy optical theorem");
        near(metal.extinction_cross_section_m2/(pi*radius*radius), metal.extinction_efficiency, 1e-15, "SI cross sections");

        rejected([&]{ (void)rayleigh_sphere(lambda, radius, 1, {1,-0.1}); }, "reject gain");
        rejected([&]{ (void)rayleigh_sphere(lambda, 0, 1, {2,0}); }, "reject zero radius");
        rejected([&]{ (void)rayleigh_sphere(lambda, 1e-6, 1, {2,0}); }, "reject large sphere");
        rejected([&]{ (void)rayleigh_sphere(lambda, radius, 1, {100,0}); }, "reject large internal size");
        rejected([&]{ (void)rayleigh_sphere(lambda, radius, 1, {-2,0}); }, "reject electrostatic pole");
        std::cout << "NOVA restricted Rayleigh sphere reference cases passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
