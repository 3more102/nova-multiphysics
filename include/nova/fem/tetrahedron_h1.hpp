#pragma once

#include "nova/fem/mesh3d.hpp"

#include <array>
#include <cstddef>

namespace nova::fem {

struct TetrahedronP1BasisEvaluation3D {
    std::array<double, 4> values{};
    std::array<std::array<double, 3>, 4> reference_gradients{};
};

// P1 nodal Lagrange basis on the reference tetrahedron with vertices
// (0,0,0), (1,0,0), (0,1,0), and (0,0,1).
[[nodiscard]] TetrahedronP1BasisEvaluation3D
evaluate_tetrahedron_p1_basis_3d(double r, double s, double t);

struct TetrahedronP1PhysicalGradients3D {
    std::array<std::array<double, 3>, 4> gradients{};
    double signed_jacobian_determinant{};
};

// Physical P1 gradients are constant over a linear tetrahedron.
[[nodiscard]] TetrahedronP1PhysicalGradients3D
evaluate_tetrahedron_p1_physical_gradients_3d(
    const TetrahedralMesh3D& mesh,
    std::size_t element
);

} // namespace nova::fem
