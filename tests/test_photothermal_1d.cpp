#include "nova/physics/photothermal_1d.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void near(double actual, double expected, double tol = 1e-10) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tol)
        throw std::runtime_error("photothermal reference mismatch");
}
}
int main() {
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
