#pragma once
#include <complex>
#include <span>
#include <vector>

namespace nova::optics {
// Electric field normalized to unit incident tangential E. SI coordinates,
// exp(-i omega t), nonmagnetic isotropic media, normal incidence.
struct Helmholtz1DResult {
    std::vector<std::complex<double>> electric_field;
    std::vector<double> cell_absorptance;
    std::complex<double> reflection;
    std::complex<double> transmission;
    double reflectance;
    double transmittance;
    double absorptance;
    double power_balance_residual;
};
// P1 FEM on an ordered mesh; one passive relative permittivity per cell.
// Exact outgoing impedance conditions in lossless exterior half spaces.
// Reference planes are the first and last mesh nodes.
[[nodiscard]] Helmholtz1DResult solve_helmholtz_1d(
    double wavelength_m, double incident_index, double exit_index,
    std::span<const double> nodes_m,
    std::span<const std::complex<double>> relative_permittivity);
} // namespace nova::optics
