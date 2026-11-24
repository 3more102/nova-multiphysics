#pragma once
#include "nova/physics/optics.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace nova::optics {
// Power-wave amplitudes at the outer faces, in left/right port order.
// The extrema are absorbed fractions for unit total incident power.
struct CoherentTwoPort {
    Complex s11, s12, s21, s22;
    double minimum_absorptance, maximum_absorptance;
};
struct CoherentPower {
    Complex outgoing_left, outgoing_right;
    double incident_power, outgoing_power, absorbed_power;
};

[[nodiscard]] inline CoherentTwoPort coherent_two_port(
    double wavelength, double left_index, double right_index,
    double left_angle, Polarization polarization, std::span<const Layer> layers) {
    if (!std::isfinite(left_index) || !std::isfinite(right_index) ||
        left_index <= 0 || right_index <= 0 || !std::isfinite(left_angle) ||
        left_angle < 0 || left_angle >= std::numbers::pi_v<double>/2)
        throw std::invalid_argument("invalid two-port exterior");
    const double ratio = (left_index / right_index) * std::sin(left_angle);
    if (!std::isfinite(ratio) || ratio >= 1)
        throw std::domain_error("two propagating exterior ports required");
    const double right_angle = std::asin(ratio);
    const auto forward = oblique_incidence_stack(wavelength, left_index, right_index,
                                                 left_angle, polarization, layers);
    const std::vector<Layer> reversed(layers.rbegin(), layers.rend());
    const auto backward = oblique_incidence_stack(wavelength, right_index, left_index,
                                                  right_angle, polarization, reversed);
    const double c1 = std::cos(left_angle), c2 = std::cos(right_angle);
    const double y1 = polarization == Polarization::te ? left_index*c1 : left_index/c1;
    const double y2 = polarization == Polarization::te ? right_index*c2 : right_index/c2;
    const double scale = std::sqrt(y2)/std::sqrt(y1);
    CoherentTwoPort s{forward.reflection, backward.transmission/scale,
                      forward.transmission*scale, backward.reflection, 0, 0};
    const double a = std::norm(s.s11)+std::norm(s.s21);
    const double d = std::norm(s.s12)+std::norm(s.s22);
    const Complex b = std::conj(s.s11)*s.s12+std::conj(s.s21)*s.s22;
    const double mid = (a+d)/2;
    const double radius = std::hypot((a-d)/2, std::abs(b));
    const double low = 1-(mid+radius), high = 1-(mid-radius);
    constexpr double tolerance = 1e-8;
    if (!std::isfinite(low) || !std::isfinite(high) ||
        std::abs(s.s12-s.s21) > tolerance || low < -tolerance || high > 1+tolerance)
        throw std::runtime_error("two-port reciprocity or passivity audit failed");
    s.minimum_absorptance = std::clamp(low,0.0,1.0);
    s.maximum_absorptance = std::clamp(high,0.0,1.0);
    return s;
}

// Incident amplitudes have units sqrt(power), or normalized equivalents.
// Absorption is a scattering residual, not a volume-integrated loss.
[[nodiscard]] inline CoherentPower coherent_power(const CoherentTwoPort& s,
                                                   Complex left, Complex right) {
    const Complex l = s.s11*left+s.s12*right;
    const Complex r = s.s21*left+s.s22*right;
    const double in = std::norm(left)+std::norm(right);
    const double out = std::norm(l)+std::norm(r);
    if (!std::isfinite(in) || !std::isfinite(out))
        throw std::overflow_error("nonfinite coherent power");
    if (out > in && out-in > 1e-8*in)
        throw std::runtime_error("coherent power violates passivity");
    return {l,r,in,out,in-out};
}
} // namespace nova::optics
