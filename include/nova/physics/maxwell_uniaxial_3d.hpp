#pragma once

#include "nova/physics/maxwell_3d.hpp"

namespace nova::physics {

// Reciprocal uniaxial permittivity: eps_o I + (eps_e-eps_o) a a^T.
// The optic axis is a finite nonzero real vector; its magnitude is ignored.
// Convention exp(-i omega t): both principal losses must be nonnegative.
// Principal real parts may be negative (driven problem only).
struct MaxwellUniaxialMaterial3D {
    std::complex<double> ordinary_permittivity{1.0,0.0};
    std::complex<double> extraordinary_permittivity{1.0,0.0};
    fem::EdgeVector3D optic_axis{0.0,0.0,1.0};
    double relative_permeability{1.0};
};

// Uses the same global edge ordering, PEC reduction and real-block
// representation as the isotropic Maxwell solver.
[[nodiscard]] MaxwellSystem3D assemble_maxwell_uniaxial_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellUniaxialMaterial3D> materials,
    double vacuum_wavenumber, bool pec_on_exterior);

// Integrate the passive dielectric loss exactly for affine Nedelec fields.
// Result is one absorbed power in watts per element.
[[nodiscard]] std::vector<double> integrate_maxwell_uniaxial_absorption_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellUniaxialMaterial3D> materials,
    double vacuum_wavenumber, const MaxwellSystem3D& system,
    std::span<const std::complex<double>> edge_field);

} // namespace nova::physics
