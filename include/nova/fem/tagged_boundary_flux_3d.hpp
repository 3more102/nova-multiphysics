#pragma once

#include "nova/fem/p1_flux_3d.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace nova::fem {

// Signed heat flow derived from each boundary owner's raw P1 element flux.
// Positive is outward. This is a diagnostic, not an equilibrated H(div) flux:
// local conservation is not guaranteed for arbitrary sources/material fields.
struct TaggedBoundaryFlux3D {
    int physical_tag{};
    std::size_t face_count{};
    double area{};               // m^2 when geometry is in metres
    double outward_heat_rate{};  // W when flux is in W/m^2
};

// The exterior-face list may be in any order or winding, but must contain only
// unique actual exterior faces. Flux entries follow mesh.elements ordering.
// Tags (including 0 for untagged faces) follow exterior_faces ordering.
// Validates topology, face ownership, field lengths, and finite arithmetic.
[[nodiscard]] std::vector<TaggedBoundaryFlux3D> evaluate_tagged_boundary_flux_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const BoundaryFace3D> exterior_faces,
    std::span<const int> physical_tags,
    std::span<const P1TetrahedralFlux3D> cell_fields
);

} // namespace nova::fem
