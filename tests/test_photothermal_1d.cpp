#include "nova/physics/photothermal_1d.hpp"
#include <algorithm>
#include <array>
#include <complex>
#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void near(double actual, double expected, double tol = 1e-10) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tol)
        throw std::runtime_error("photothermal reference mismatch");
}

void end_to_end_helmholtz() {
    using namespace nova::optics;
    constexpr std::size_t cells = 120;
    constexpr double length = 1.2e-6, conductivity = 0.5, irradiance = 2e7;
    constexpr double left = 300, right = 304;
    std::vector<double> x(cells + 1), k(cells, conductivity);
    std::vector<std::complex<double>> eps(cells);
    for (std::size_t j = 0; j <= cells; ++j) {
        const double u = static_cast<double>(j) / static_cast<double>(cells);
        x[j] = length * (0.5 * u + 0.5 * u * u);
        if (j < cells) eps[j] = j < cells / 2
            ? std::complex<double>{2.25, 0.0} : std::complex<double>{4.0, 0.8};
    }
    const auto optical = solve_helmholtz_1d(1.55e-6, 1.0, 1.4, x, eps);
    if (!(optical.absorptance > 0.01 && optical.absorptance < 1.0))
        throw std::runtime_error("expected finite nontrivial optical absorption");
    near(optical.power_balance_residual, 0, 1e-9);
    const auto thermal = solve_photothermal_1d(x, k, optical, irradiance, left, right);
    near(thermal.absorbed_power_W_m2 / irradiance, optical.absorptance, 1e-12);
    near(thermal.energy_residual_W_m2 / irradiance, 0, 1e-12);
    double integrated = 0;
    for (std::size_t j = 0; j < cells; ++j) {
        integrated += thermal.cell_heat_W_m3[j] * (x[j+1] - x[j]);
        if (j < cells / 2) near(thermal.cell_heat_W_m3[j], 0);
    }
    near(integrated / irradiance, optical.absorptance, 1e-12);
    // Independent fixed-temperature Green function for -k*T''=q.
    // At each mesh node the kernel is affine across each source cell,
    // so midpoint integration of the cell-constant heating is exact.
    for (std::size_t i = 0; i <= cells; ++i) {
        double expected = left + (right - left) * x[i] / length;
        for (std::size_t j = 0; j < cells; ++j) {
            const double s = 0.5 * (x[j] + x[j+1]);
            const double g = std::min(x[i], s) *
                (length - std::max(x[i], s)) / length;
            expected += g * irradiance * optical.cell_absorptance[j] / conductivity;
        }
        near(thermal.temperature_K[i], expected, 1e-8);
    }
    const auto dark = solve_photothermal_1d(x, k, optical, 0, left, right);
    const auto doubled = solve_photothermal_1d(x, k, optical, 2 * irradiance, left, right);
    for (std::size_t i = 0; i <= cells; ++i) {
        const double baseline = left + (right - left) * x[i] / length;
        near(dark.temperature_K[i], baseline, 1e-8);
        near(doubled.temperature_K[i] - baseline,
            2 * (thermal.temperature_K[i] - baseline), 1e-8);
    }
    // A genuinely lossless optical solve must yield no deposited heat.
    for (auto& material : eps) material = {material.real(), 0.0};
    const auto transparent = solve_helmholtz_1d(1.55e-6, 1.0, 1.4, x, eps);
    const auto unheated = solve_photothermal_1d(x, k, transparent, irradiance, left, right);
    near(unheated.absorbed_power_W_m2, 0);
    for (std::size_t i = 0; i <= cells; ++i)
        near(unheated.temperature_K[i], dark.temperature_K[i], 1e-8);
}
}
int main() {
    end_to_end_helmholtz();
    using nova::optics::solve_photothermal_1d;
    const std::array<double, 3> x{0, 0.5, 1};
    const std::array<double, 2> k{2, 2}, a{0.25, 0.25}, zero{0, 0};
    const auto symmetric = solve_photothermal_1d(x, k, a, 80, 300, 300);
    near(symmetric.temperature_K[1], 302.5);
    near(symmetric.temperature_K[2], 300);
    near(symmetric.face_flux_W_m2[0], -20);
    near(symmetric.face_flux_W_m2[1], 0);
    near(symmetric.face_flux_W_m2[2], 20);
    near(symmetric.cell_heat_W_m3[0], 40);
    near(symmetric.energy_residual_W_m2, 0);
    const std::array<double, 2> unequal_k{1, 3};
    const auto conduction = solve_photothermal_1d(x, unequal_k, zero, 80, 320, 300);
    near(conduction.face_flux_W_m2[0], 30);
    near(conduction.temperature_K[1], 305);
    near(conduction.temperature_K[2], 300);
    const std::array<double, 2> one_sided{0.5, 0};
    const auto asymmetric = solve_photothermal_1d(x, unequal_k, one_sided, 80, 300, 300);
    near(asymmetric.face_flux_W_m2[0], -25);
    near(asymmetric.face_flux_W_m2[2], 15);
    near(asymmetric.temperature_K[1], 302.5);
    near(asymmetric.temperature_K[2], 300);
    nova::optics::Helmholtz1DResult optical{};
    optical.cell_absorptance = {0.25, 0.25};
    const auto coupled = solve_photothermal_1d(x, k, optical, 80, 300, 300);
    near(coupled.temperature_K[1], 302.5);
    const std::array<double, 2> invalid_a{0.8, 0.8};
    bool rejected = false;
    try { (void)solve_photothermal_1d(x, k, invalid_a, 80, 300, 300); }
    catch (const std::invalid_argument&) { rejected = true; }
    if (!rejected) return 1;
    std::cout << "photothermal 1D references passed\n";
}
