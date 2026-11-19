#pragma once

#include "nova/fem/diffusion_3d.hpp"

#include <span>
#include <vector>

namespace nova::fem {

// Per-cell physical gradient and element-averaged diffusive flux q = -K grad(u).
// Entries follow TetrahedralMesh3D::elements ordering. P1 gradients are constant
// per cell; the material is sampled using the same positive degree-two
// tetrahedral quadrature rule as the NOVA 3D FEM assembly.
struct P1TetrahedralFlux3D {
    Point3D gradient;
    Point3D flux;
};

// No solver is run. The caller supplies its P1 nodal field (e.g. a stationary
// diffusion result). Uses the same conductivity selection as the solve.
// Throws on invalid topology, nonfinite nodal data or invalid conductivity.
[[nodiscard]] std::vector<P1TetrahedralFlux3D> evaluate_p1_tetrahedral_flux_3d(
    const ScalarDiffusion3DProblem& problem,
    std::span<const double> nodal_values
);

} // namespace nova::fem
