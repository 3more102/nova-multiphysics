#pragma once
#include <complex>
#include <span>

namespace nova::optics {
using Complex = std::complex<double>;
// SI units; exp(-i omega t), passive media have Im(epsilon), Im(n) >= 0.
struct Layer {
    Complex index;
    double thickness_m;
};
struct StackResult {
    Complex reflection;
    Complex transmission;
    double reflectance;
    double transmittance;
    double absorptance;
};
// Normal incidence, nonmagnetic isotropic layers, real positive exterior indices.
[[nodiscard]] StackResult normal_incidence_stack(double wavelength_m,
    double incident_index, double exit_index, std::span<const Layer> layers);
// Passive Drude relative permittivity; frequencies in rad/s.
[[nodiscard]] Complex drude_permittivity(double omega, double epsilon_infinity,
    double plasma_frequency, double collision_frequency);
[[nodiscard]] Complex passive_index(Complex relative_permittivity);

struct KerrResult {
    double index_shift;
    double phase_shift_rad;
};
// Prescribed uniform intensity, instantaneous scalar Kerr response; not a self-consistent solve.
[[nodiscard]] KerrResult kerr_phase(double wavelength_m, double n2_m2_per_w,
    double intensity_w_per_m2, double length_m);

struct SurfacePlasmon {
    Complex beta_per_m;
    Complex dielectric_decay_per_m;
    Complex metal_decay_per_m;
    double intensity_propagation_length_m;
};
// Bound TM mode at a planar interface between a lossless dielectric and passive metal.
[[nodiscard]] SurfacePlasmon planar_surface_plasmon(double wavelength_m,
    double dielectric_permittivity, Complex metal_permittivity);
} // namespace nova::optics
