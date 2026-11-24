#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace nova::optics {

// Uniform, scalar, continuous-wave intensity model in SI units.
// dI/dz = -alpha*I - beta*I^2; dphi/dz = (2*pi/lambda)*n2*I.
// No diffraction, dispersion, gain, carrier dynamics or field redistribution.
struct NonlinearPropagationResult {
    double output_intensity_w_per_m2;
    double integrated_intensity_w_per_m;
    double kerr_phase_rad;
    double absorbed_intensity_w_per_m2;
};

[[nodiscard]] inline NonlinearPropagationResult nonlinear_propagation(
    double wavelength_m, double input_intensity_w_per_m2,
    double length_m, double linear_loss_per_m,
    double two_photon_absorption_m_per_w, double kerr_n2_m2_per_w) {
    const auto nonnegative = [](double x) {
        return std::isfinite(x) && x >= 0.0;
    };
    if (!std::isfinite(wavelength_m) || wavelength_m <= 0.0 ||
        !nonnegative(input_intensity_w_per_m2) || !nonnegative(length_m) ||
        !nonnegative(linear_loss_per_m) ||
        !nonnegative(two_photon_absorption_m_per_w) ||
        !std::isfinite(kerr_n2_m2_per_w))
        throw std::invalid_argument("invalid nonlinear propagation input");
    if (length_m == 0.0 || input_intensity_w_per_m2 == 0.0)
        return {input_intensity_w_per_m2, 0.0, 0.0, 0.0};

    const double optical_depth = linear_loss_per_m * length_m;
    if (!std::isfinite(optical_depth))
        throw std::overflow_error("nonfinite optical depth");
    const double effective_length = optical_depth == 0.0 ? length_m :
        -std::expm1(-optical_depth) / linear_loss_per_m;
    const double fluence_z = input_intensity_w_per_m2 * effective_length;
    const double q = two_photon_absorption_m_per_w * fluence_z;
    if (!std::isfinite(fluence_z) || !std::isfinite(q))
        throw std::overflow_error("nonfinite nonlinear optical depth");

    const double integral = fluence_z * (q == 0.0 ? 1.0 : std::log1p(q) / q);
    const double output = input_intensity_w_per_m2 * std::exp(-optical_depth) / (1.0 + q);
    // Evaluate the absorbed fraction directly to avoid cancellation for weak loss.
    const double absorbed = input_intensity_w_per_m2 *
        ((-std::expm1(-optical_depth) + q) / (1.0 + q));
    const double phase = kerr_n2_m2_per_w == 0.0 ? 0.0 :
        (2.0 * std::numbers::pi_v<double>) *
        (kerr_n2_m2_per_w * integral / wavelength_m);
    if (!std::isfinite(integral) || !std::isfinite(output) ||
        !std::isfinite(absorbed) || !std::isfinite(phase))
        throw std::overflow_error("nonfinite nonlinear propagation result");
    return {output, integral, phase, absorbed};
}

} // namespace nova::optics
