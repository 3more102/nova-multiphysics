#pragma once

#include <cmath>
#include <complex>
#include <stdexcept>

namespace nova::optics {

struct DirectionalCouplerResult {
    std::complex<double> output_a;
    std::complex<double> output_b;
    double input_power_w;
    double output_power_w;
    double dissipated_power_w;
};

// Uniform reciprocal two-mode coupled-mode reference. Amplitudes have units
// sqrt(W); kappa and delta_beta are real in rad/m, alpha is a common power
// attenuation coefficient in 1/m. The average propagation phase is removed.
// exp(-i omega t); dA/dz = i H A - alpha A/2,
// H = [[delta_beta/2, kappa], [kappa, -delta_beta/2]].
[[nodiscard]] inline DirectionalCouplerResult directional_coupler(
    std::complex<double> input_a, std::complex<double> input_b,
    double kappa, double delta_beta, double alpha, double length_m) {
    const auto finite = [](std::complex<double> z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    };
    if (!finite(input_a) || !finite(input_b) || !std::isfinite(kappa) ||
        !std::isfinite(delta_beta) || !std::isfinite(alpha) ||
        !std::isfinite(length_m) || alpha < 0.0 || length_m < 0.0)
        throw std::invalid_argument("invalid directional coupler parameters");
    const double pin = std::norm(input_a) + std::norm(input_b);
    if (!std::isfinite(pin))
        throw std::overflow_error("directional coupler input power overflow");
    if (length_m == 0.0)
        return {input_a, input_b, pin, pin, 0.0};
    const double d = delta_beta / 2.0;
    const double g = std::hypot(kappa, d);
    const double phase = g * length_m;
    if (!std::isfinite(phase))
        throw std::overflow_error("directional coupler phase overflow");
    const double loss = alpha * length_m;
    const double attenuation = std::exp(-0.5 * loss);
    const double c = std::cos(phase);
    const double s = std::sin(phase);
    const double u = g == 0.0 ? 0.0 : s * (d / g);
    const double v = g == 0.0 ? 0.0 : s * (kappa / g);
    const auto a = attenuation * (std::complex<double>{c, u} * input_a +
                                    std::complex<double>{0.0, v} * input_b);
    const auto b = attenuation * (std::complex<double>{0.0, v} * input_a +
                                    std::complex<double>{c, -u} * input_b);
    const double pout = std::norm(a) + std::norm(b);
    if (!finite(a) || !finite(b) || !std::isfinite(pout))
        throw std::overflow_error("directional coupler output overflow");
    return {a, b, pin, pout, -std::expm1(-loss) * pin};
}

} // namespace nova::optics
