#pragma once

#include "nova/physics/electrostatics_3d.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace nova::physics {

// Approximate charge on a prescribed-potential conductor, using the P1
// element-averaged electric displacement on the adjacent dielectric side.
// Positive conductor charge points into the dielectric: Q = -integral(D.n) dA,
// where n is the outward normal of the DIELECTRIC computational domain.
// This raw boundary-flux diagnostic is not locally equilibrated or guaranteed
// to obey Gauss's law on coarse/unstructured meshes.
struct ElectrodeCharge3D {
    int electrode_id{};
    std::size_t face_count{};
    double area_square_metres{};
    double potential_volts{};
    double dielectric_outward_charge_coulombs{};
    double conductor_charge_coulombs{};
};

// Faces must be unique exterior triangles with vertices all fixed by the
// problem's prescribed-potential conditions. Every face belonging to one
// electrode ID must have a common potential to floating-point precision.
// The returned groups are sorted by nonnegative electrode ID. Winding is
// ignored; geometric normals are oriented outward from the dielectric.
// Neither the fields nor the mesh are re-solved or modified.
[[nodiscard]] std::vector<ElectrodeCharge3D> evaluate_electrode_charges_3d(
    const Electrostatics3DProblem& problem,
    const Electrostatics3DResult& solution,
    std::span<const fem::BoundaryFace3D> electrode_faces,
    std::span<const int> electrode_ids
);

} // namespace nova::physics
