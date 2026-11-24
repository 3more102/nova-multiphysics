#pragma once
#include <cstddef>
#include <vector>

namespace nova::physics {

enum class SlabPolarization { te, tm };

// Symmetric, lossless, isotropic, nonmagnetic planar dielectric slab.
// All lengths are in metres. The core occupies |x| < half_thickness.
struct SlabWaveguide {
    double vacuum_wavelength;
    double half_thickness;
    double core_index;
    double cladding_index;
};

struct SlabGuidedMode {
    std::size_t order; // 0 is fundamental; even orders have even scalar field
    SlabPolarization polarization;
    double effective_index;
    double propagation_constant; // rad/m
    double cladding_decay; // 1/m
    double normalized_transverse; // u
    double normalized_decay; // w
    double core_power_fraction; // integrated longitudinal Poynting flux
    double dispersion_residual; // radians
};

// Exact slab dispersion reference, not a discretized FEM eigenmode solver.
// Only strictly guided modes are returned; cutoff modes are not bound.
// Rejects unsupported numerical scales and unrepresentable near-cutoff modes.
[[nodiscard]] std::vector<SlabGuidedMode> solve_symmetric_slab(
    const SlabWaveguide& guide, SlabPolarization polarization);

} // namespace nova::physics
