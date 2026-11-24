#pragma once
#include "nova/physics/optics.hpp"
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

// Scalar field normalized to unit incident scalar amplitude: E_y for TE,
// Z0*H_y for TM (the TM scalar is normalized independently of E_x).
// Reflection/transmission use tangential E amplitudes, matching optics.hpp.
struct ObliqueHelmholtz1DResult {
    std::vector<std::complex<double>> scalar_field;
    std::vector<double> cell_absorptance;
    std::complex<double> reflection;
    std::complex<double> transmission;
    double reflectance;
    double transmittance;
    double absorptance;
    double power_balance_residual;
};
// Isotropic nonmagnetic passive layers; real positive exterior indices.
// Angle from normal in radians. Evanescent and exact critical-angle exits
// are supported (zero normal transmitted power); zero TM permittivity excluded.
[[nodiscard]] ObliqueHelmholtz1DResult solve_oblique_helmholtz_1d(
    double wavelength_m, double incident_index, double exit_index,
    double angle_rad, Polarization polarization, std::span<const double> nodes_m,
    std::span<const std::complex<double>> relative_permittivity);
} // namespace nova::optics
