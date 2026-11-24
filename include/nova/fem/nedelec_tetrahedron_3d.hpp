#pragma once

#include "nova/fem/mesh3d.hpp"

#include <array>

namespace nova::fem {

using Edge3D = std::array<std::size_t, 2>;
using EdgeVector3D = std::array<double, 3>;
using EdgeMatrix6 = std::array<std::array<double, 6>, 6>;

inline constexpr std::array<Edge3D, 6> tetrahedron_local_edges{{
    {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}
}};

struct NedelecTetrahedron3D {
    std::array<Edge3D, 6> edges{}; // globally ascending node indices
    std::array<EdgeVector3D, 4> gradients{};
    std::array<double, 6> signs{};
    std::array<EdgeVector3D, 6> curls{};
    double volume{};
};

// Lowest-order first-family H(curl) basis, normalized to unit oriented
// edge circulation. Supports either tetrahedron orientation.
[[nodiscard]] NedelecTetrahedron3D make_nedelec_tetrahedron_3d(
    const TetrahedralMesh3D& mesh, std::size_t element);

[[nodiscard]] std::array<EdgeVector3D, 6> evaluate_nedelec_basis_3d(
    const NedelecTetrahedron3D& element, double r, double s, double t);

struct NedelecElementMatrices3D {
    EdgeMatrix6 curl_curl{};
    EdgeMatrix6 mass{};
};

// Exact integrals for constant isotropic coefficients on an affine cell.
[[nodiscard]] NedelecElementMatrices3D nedelec_element_matrices_3d(
    const NedelecTetrahedron3D& element);

} // namespace nova::fem
