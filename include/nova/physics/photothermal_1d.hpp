#pragma once
#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>
#include "nova/physics/optical_helmholtz_1d.hpp"

namespace nova::optics {
// Steady planar conduction with fixed face temperatures, piecewise-constant
// conductivity and uniform volumetric heat within each optical cell.
// Flux is positive in +x; all fluxes are per unit transverse area.
struct Photothermal1DResult {
    std::vector<double> temperature_K;
    std::vector<double> face_flux_W_m2;
    std::vector<double> cell_heat_W_m3;
    double absorbed_power_W_m2{};
    double energy_residual_W_m2{};
};

[[nodiscard]] inline Photothermal1DResult solve_photothermal_1d(
    std::span<const double> nodes_m,
    std::span<const double> conductivity_W_mK,
    std::span<const double> cell_absorptance,
    double incident_flux_W_m2, double left_K, double right_K) {
    const auto n = conductivity_W_mK.size();
    if (n == 0 || nodes_m.size() != n + 1 || cell_absorptance.size() != n ||
        !std::isfinite(incident_flux_W_m2) || incident_flux_W_m2 < 0 ||
        !std::isfinite(left_K) || !std::isfinite(right_K) || left_K < 0 || right_K < 0)
        throw std::invalid_argument("invalid photothermal inputs");
    Photothermal1DResult out;
    out.temperature_K.resize(n + 1);
    out.face_flux_W_m2.resize(n + 1);
    out.cell_heat_W_m3.resize(n);
    std::vector<double> resistance(n), power(n);
    double total_r = 0, weighted_heat = 0, accumulated = 0, total_a = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double h = nodes_m[i + 1] - nodes_m[i];
        const double k = conductivity_W_mK[i], a = cell_absorptance[i];
        if (!std::isfinite(h) || h <= 0 || !std::isfinite(k) || k <= 0 ||
            !std::isfinite(a) || a < 0 || a > 1)
            throw std::invalid_argument("invalid optical/thermal cell");
        resistance[i] = h / k;
        power[i] = incident_flux_W_m2 * a;
        out.cell_heat_W_m3[i] = power[i] / h;
        weighted_heat += resistance[i] * (accumulated + 0.5 * power[i]);
        total_r += resistance[i];
        accumulated += power[i];
        total_a += a;
        if (!std::isfinite(out.cell_heat_W_m3[i]) ||
            !std::isfinite(weighted_heat) || !std::isfinite(total_r) ||
            !std::isfinite(accumulated) || resistance[i] <= 0)
            throw std::overflow_error("photothermal arithmetic out of range");
    }
    if (total_a > 1.0 + 1e-8)
        throw std::invalid_argument("cell absorptance exceeds incident power");
    out.absorbed_power_W_m2 = accumulated;
    out.temperature_K[0] = left_K;
    out.face_flux_W_m2[0] = (left_K - right_K - weighted_heat) / total_r;
    for (std::size_t i = 0; i < n; ++i) {
        out.temperature_K[i + 1] = out.temperature_K[i] -
            resistance[i] * (out.face_flux_W_m2[i] + 0.5 * power[i]);
        out.face_flux_W_m2[i + 1] = out.face_flux_W_m2[i] + power[i];
        if (!std::isfinite(out.temperature_K[i + 1]) ||
            !std::isfinite(out.face_flux_W_m2[i + 1]))
            throw std::overflow_error("photothermal solution out of range");
    }
    out.energy_residual_W_m2 = out.face_flux_W_m2.back() -
        out.face_flux_W_m2.front() - accumulated;
    if (!std::isfinite(out.energy_residual_W_m2))
        throw std::overflow_error("photothermal energy residual out of range");
    return out;
}

[[nodiscard]] inline Photothermal1DResult solve_photothermal_1d(
    std::span<const double> nodes_m,
    std::span<const double> conductivity_W_mK,
    const Helmholtz1DResult& optical,
    double incident_flux_W_m2, double left_K, double right_K) {
    return solve_photothermal_1d(nodes_m, conductivity_W_mK,
        optical.cell_absorptance, incident_flux_W_m2, left_K, right_K);
}
} // namespace nova::optics
