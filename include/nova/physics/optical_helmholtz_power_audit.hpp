#pragma once
#include "nova/physics/optical_helmholtz_1d.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace nova::optics {
struct HelmholtzCellPower {
    double left_flux;
    double right_flux;
    double absorption_residual;
};
struct HelmholtzPowerAudit {
    std::vector<HelmholtzCellPower> cells;
    double max_interface_jump;
    double entrance_residual;
    double exit_residual;
};

// Signed rightward flux, normalized to incident power. Uses discrete
// element endpoint reactions, not raw piecewise-constant field gradients.
// Inputs must describe the same normal-incidence solve.
[[nodiscard]] inline HelmholtzPowerAudit audit_helmholtz_cell_power(
    double wavelength_m, double incident_index,
    std::span<const double> nodes_m,
    std::span<const std::complex<double>> relative_permittivity,
    const Helmholtz1DResult& solution) {
    const auto n = relative_permittivity.size();
    if (!std::isfinite(wavelength_m) || wavelength_m <= 0 ||
        !std::isfinite(incident_index) || incident_index <= 0 ||
        n == 0 || nodes_m.size() - (nodes_m.empty() ? 0u : 1u) != n ||
        solution.electric_field.size() != nodes_m.size() ||
        solution.cell_absorptance.size() != n ||
        !std::isfinite(solution.reflectance) || !std::isfinite(solution.transmittance))
        throw std::invalid_argument("invalid Helmholtz power audit inputs");
    const double k = 2.0 * std::numbers::pi_v<double> / wavelength_m;
    HelmholtzPowerAudit audit{};
    audit.cells.reserve(n);
    for (std::size_t j = 0; j < n; ++j) {
        const double h = k * (nodes_m[j+1] - nodes_m[j]);
        const auto e = relative_permittivity[j];
        const auto a = solution.electric_field[j];
        const auto b = solution.electric_field[j+1];
        const double loss = solution.cell_absorptance[j];
        if (!std::isfinite(h) || h <= 0 ||
            !std::isfinite(e.real()) || !std::isfinite(e.imag()) || e.imag() < 0 ||
            !std::isfinite(loss) || loss < 0)
            throw std::invalid_argument("invalid Helmholtz audit cell");
        const auto ab = std::conj(a)*b;
        const double gradient_flux = std::imag(ab)/h;
        const double left = (gradient_flux + h*std::imag(e*(std::norm(a)/3.0 + ab/6.0)))/incident_index;
        const double right = (gradient_flux - h*std::imag(e*(std::norm(b)/3.0 + std::conj(ab)/6.0)))/incident_index;
        const double residual = left-right-loss;
        if (!std::isfinite(left) || !std::isfinite(right) || !std::isfinite(residual))
            throw std::overflow_error("nonfinite Helmholtz audit flux");
        if (j != 0)
            audit.max_interface_jump = std::max(audit.max_interface_jump,
                std::abs(left-audit.cells.back().right_flux));
        audit.cells.push_back({left, right, residual});
    }
    audit.entrance_residual = audit.cells.front().left_flux-(1.0-solution.reflectance);
    audit.exit_residual = audit.cells.back().right_flux-solution.transmittance;
    if (!std::isfinite(audit.entrance_residual) || !std::isfinite(audit.exit_residual))
        throw std::overflow_error("nonfinite Helmholtz audit endpoint");
    return audit;
}
} // namespace nova::optics
