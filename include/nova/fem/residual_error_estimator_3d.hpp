#pragma once

#include "nova/fem/interior_flux_jumps_3d.hpp"

#include <functional>
#include <span>
#include <vector>

namespace nova::fem {

// A selected natural-boundary face with prescribed outward conductive flux.
// Include zero-flux exterior faces explicitly when their residual is wanted.
// Dirichlet faces must not be included.
struct NeumannResidualFace3D {
    BoundaryFace3D face{};
    std::function<double(double, double, double)> outward_flux;
};

// Robin natural face: q.n = h(x) * (u_h(x) - ambient(x)).
struct RobinResidualFace3D {
    BoundaryFace3D face{};
    std::function<double(double, double, double)> transfer_coefficient;
    std::function<double(double, double, double)> ambient_value;
};

enum class BoundaryResidualKind3D { Neumann, Robin };

struct BoundaryFluxResidualFace3D {
    BoundaryFace3D nodes{}; // sorted exterior vertex IDs
    std::size_t owner_cell{};
    double area{};
    double diameter{};
    double mismatch_l2{}; // L2 mismatch to the selected natural BC
    BoundaryResidualKind3D kind{BoundaryResidualKind3D::Neumann};
};

struct DiffusionResidualIndicators3D {
    // Cell-source residual contribution h_K * ||f||_{L2(K)}.
    std::vector<double> cell_residual_errors;
    // Existing interior normal-flux jump report, shared equally by adjacent cells.
    InteriorFluxJumpReport3D flux_jumps;
    // Nonnegative Neumann mismatch contribution per tetrahedron.
    std::vector<double> boundary_residual_errors;
    // Only supplied exterior natural faces are included in this report.
    std::vector<BoundaryFluxResidualFace3D> boundary_residual_faces;
    // Combined nonnegative element indicators: sqrt(cell^2 + face_jump^2 + boundary^2).
    std::vector<double> element_errors;
    // Euclidean aggregate: sqrt(sum_K eta_K^2).
    double total_error{};
};

// P1 residual indicator for -div(k grad(u_h)) = f on tetrahedra.
// Conductivity is positive and constant within each element. It combines the
// cell-source residual with interior normal-flux jumps. Natural-boundary
// residuals, within-cell coefficient gradients and tensor materials are
// omitted by this four-argument overload; neither overload is a certified
// error bound. No unlabeled or Dirichlet boundary face is automatically added.
[[nodiscard]] DiffusionResidualIndicators3D estimate_diffusion_residual_error_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const double> nodal_values,
    std::span<const double> element_conductivity,
    const std::function<double(double, double, double)>& source
);

// As above, but also include h_F * ||q_h dot n - g_N||^2_L2(F) for
// explicitly selected natural exterior faces, once per owning tetrahedron.
// Face winding is irrelevant; duplicate, interior or invalid faces are
// rejected. A positive six-point degree-four triangle quadrature integrates
// the prescribed boundary mismatch. Unlisted faces contribute no boundary
// residual; include insulated faces with a zero callback as needed.
[[nodiscard]] DiffusionResidualIndicators3D estimate_diffusion_residual_error_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const double> nodal_values,
    std::span<const double> element_conductivity,
    const std::function<double(double, double, double)>& source,
    std::span<const NeumannResidualFace3D> natural_boundary_faces
);

// Includes disjoint explicit Neumann and Robin exterior faces; Robin flux is
// evaluated from the P1 trace at each face quadrature point. Omitted exterior
// faces are NOT automatically classified as insulated.
[[nodiscard]] DiffusionResidualIndicators3D estimate_diffusion_residual_error_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const double> nodal_values,
    std::span<const double> element_conductivity,
    const std::function<double(double, double, double)>& source,
    std::span<const NeumannResidualFace3D> neumann_boundary_faces,
    std::span<const RobinResidualFace3D> robin_boundary_faces
);

} // namespace nova::fem
