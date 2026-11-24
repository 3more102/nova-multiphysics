# Passive dispersive optical materials

The header `nova/physics/dispersive_optical_material.hpp` evaluates isotropic
relative electric permittivity at a specified real, positive angular frequency.
It uses the same exp(-i omega t) convention as the driven Maxwell FEM API.

    epsilon_r(omega) = epsilon_inf + sum_j S_j^2 /
                        (Omega_j^2 - omega^2 - i gamma_j omega)

All frequency parameters (omega, Omega, S, gamma) are in rad/s. A zero
resonance is a Drude term; a positive resonance is a Lorentz term.
epsilon_inf must be positive and real; S and gamma must be nonnegative.
These restrictions describe a passive model (Im epsilon_r >= 0).
Negative real permittivity is permitted. Zero-strength poles are skipped.
An undamped pole evaluated exactly at resonance is rejected.

## Integration with the driven Maxwell API

~~~cpp
#include <array>
#include <vector>
#include "nova/physics/dispersive_optical_material.hpp"
#include "nova/physics/maxwell_3d.hpp"

// Illustrative synthetic parameters, NOT a fitted gold or silver model.
const double omega = nova::physics::optical_angular_frequency(1.55e-6);
const std::array<nova::physics::OpticalPole, 1> poles{
    nova::physics::OpticalPole{0.0, 2.0 * omega, 0.05 * omega}};
const auto eps = nova::physics::dispersive_permittivity(omega, 1.0, poles);
const auto n = nova::physics::passive_refractive_index(eps);
// Given an existing validated tetrahedral mesh:
const std::vector<nova::physics::MaxwellCellMaterial3D> materials(
    mesh.elements.size(), {eps, 1.0});
const auto system = nova::physics::assemble_maxwell_3d(
    mesh, materials, omega / 299792458.0, false);
~~~

For heterogeneous systems evaluate each region's own model at the same
frequency and populate the material array in mesh cell order.
Reassemble for each frequency. The index helper is intended only for
nonmagnetic isotropic media and uses the passive principal square root.

## Verification gates

The added cases in test_maxwell_3d.cpp check analytical Drude and Lorentz
values, superposition, passive index, wavelength conversion, large-scale
frequency normalization, invalid damping and undamped singularities.
A manufactured negative-real-permittivity Maxwell case checks complex
edge circulations and positive integrated dielectric absorption against
the analytical constant-field value. Run the existing nova_tests target
through CTest. Test registration is unchanged.

## Scope

This is a local, scalar frequency-domain constitutive model, not a
validated metal database or a surface-plasmon propagation benchmark.
No nonlinear susceptibility, anisotropy, nonlocal response, gain,
PML, port excitation, or dispersive eigenfrequency solver is introduced.
Do not feed complex or negative epsilon into the existing lossless
positive-permittivity cavity eigensolver. Material fits require independent
validation over the intended frequency and temperature range.
