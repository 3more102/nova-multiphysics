#pragma once

#include "nova/physics/maxwell_cavity_3d.hpp"

namespace nova::physics {

struct MaxwellCavityCellEnergy3D {
    double electric_integral{}; // integral epsilon_r |E|^2 dV
    double magnetic_integral{}; // integral mu_r^-1 |curl E|^2 / k0^2 dV
    double electric_participation{};
    double total_participation{};
};

struct MaxwellCavityEnergy3D {
    std::vector<MaxwellCavityCellEnergy3D> cells;
    double electric_integral{};
    double magnetic_integral{};
    double relative_energy_imbalance{};
    double peak_electric_weight{};
    double electric_mode_volume_m3{};
};

// Lossless, nondispersive, positive scalar materials only. Edge amplitudes
// are those of the supplied real mode (arbitrary normalization). Integrals
// become time-averaged SI energies after multiplication by epsilon_0/4
// when the mesh is in metres and the electric field in V/m.
// Mode volume is integral epsilon_r |E|^2 / max epsilon_r |E|^2.
// The peak is the exact cellwise affine-field maximum, including both
// material traces at interfaces. No claim of dispersive mode normalization.
[[nodiscard]] MaxwellCavityEnergy3D evaluate_maxwell_cavity_energy_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials,
    const MaxwellCavityMatrices3D& matrices,
    const MaxwellCavityMode3D& mode);

} // namespace nova::physics
