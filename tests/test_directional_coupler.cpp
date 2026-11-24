#include "nova/physics/directional_coupler.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace {
void check(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
void near(std::complex<double> a, std::complex<double> b, const char* message) {
    check(std::abs(a-b) <= 2e-12 * std::max(1.0, std::abs(b)), message);
}
}

int main() {
    using nova::optics::directional_coupler;
    using C = std::complex<double>;
    const double pi = std::numbers::pi_v<double>;
    try {
        const auto cross = directional_coupler(1., 0., 2., 0., 0., pi/4);
        near(cross.output_a, 0., "complete transfer: bar");
        near(cross.output_b, C{0., 1.}, "complete transfer: cross");
        const auto split = directional_coupler(1., 0., 2., 0., 0., pi/8);
        near(std::norm(split.output_a), 0.5, "3 dB bar");
        near(std::norm(split.output_b), 0.5, "3 dB cross");
        const auto mismatch = directional_coupler(1., 0., 3., 8., 0., pi/10);
        near(std::norm(mismatch.output_b), 0.36, "detuned transfer ceiling");
        near(mismatch.output_power_w, 1., "detuned conservation");
        const C a{0.2, -0.7}, b{0.9, 0.3};
        const auto zero = directional_coupler(a,b,3.,8.,2.,0.);
        near(zero.output_a,a,"zero length A");
        near(zero.output_b,b,"zero length B");
        const auto uncoupled = directional_coupler(a,b,0.,8.,0.,0.3);
        near(uncoupled.output_a, a*std::polar(1.,1.2),"uncoupled A phase");
        near(uncoupled.output_b, b*std::polar(1.,-1.2),"uncoupled B phase");
        const auto common = directional_coupler(a,b,0.,0.,2.,0.3);
        near(common.output_a,a*std::exp(-0.3),"zero Hamiltonian");
        const auto full = directional_coupler(a,b,-3.,8.,0.7,0.6);
        near(full.output_power_w,full.input_power_w*std::exp(-0.42),"loss law");
        near(full.output_power_w+full.dissipated_power_w,full.input_power_w,
             "coherent two-input energy budget");
        const auto half = directional_coupler(a,b,-3.,8.,0.7,0.3);
        const auto twice = directional_coupler(half.output_a,half.output_b,-3.,8.,0.7,0.3);
        near(twice.output_a,full.output_a,"composition A");
        near(twice.output_b,full.output_b,"composition B");
        const auto reverse = directional_coupler(0.,1.,3.,8.,0.,pi/10);
        near(reverse.output_a,mismatch.output_b,"reciprocity");
        const auto opaque = directional_coupler(a,b,3.,8.,1e300,1e10);
        near(opaque.output_power_w,0.,"large loss");
        near(opaque.dissipated_power_w,opaque.input_power_w,"large loss balance");
        const auto invalid = [](double k, double d, double alpha, double l) {
            try { (void)directional_coupler(1.,0.,k,d,alpha,l); }
            catch (const std::invalid_argument&) { return true; }
            return false;
        };
        check(invalid(1.,0.,-1.,1.),"reject gain");
        check(invalid(1.,0.,0.,-1.),"reject negative length");
        check(invalid(std::numeric_limits<double>::quiet_NaN(),0.,0.,1.),"reject NaN");
        check(invalid(1.,std::numeric_limits<double>::infinity(),0.,1.),"reject infinity");
        bool overflow = false;
        try { (void)directional_coupler(1.,0.,1e300,0.,0.,1e100); }
        catch (const std::overflow_error&) { overflow = true; }
        check(overflow,"reject phase overflow");
        std::cout << "Directional coupler reference checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
