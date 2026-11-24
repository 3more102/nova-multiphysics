#pragma once

#include "nova/fem/conductivity_3d.hpp"
#include "nova/fem/h1_space_3d.hpp"
#include "nova/linalg/conjugate_gradient.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

using ScalarField3D = SpatialScalarField3D;

struct DirichletConstraint3D {
    std::size_t node{};
    double value{};
};

struct DirichletBoundaryFace3D {
    BoundaryFace3D face{};
    ScalarField3D value;
};

struct NeumannBoundaryFace3D {
    BoundaryFace3D face{};
    // Outward physical flux q_n = -(K grad(u)) dot n.
    ScalarField3D outward_flux;
};

struct RobinBoundaryFace3D {
    BoundaryFace3D face{};
    // Outward flux q_n = h * (u - u_infinity), with finite h >= 0.
    ScalarField3D transfer_coefficient;
    ScalarField3D ambient_value;
};

struct ScalarDiffusion3DProblem {
    TetrahedralMesh3D mesh;
    ScalarField3D diffusivity;
    // Optional spatial volumetric source density; an empty callback means zero.
    ScalarField3D source;
    std::vector<DirichletConstraint3D> dirichlet_constraints;
    std::vector<DirichletBoundaryFace3D> dirichlet_faces;
    std::vector<NeumannBoundaryFace3D> neumann_faces;
    std::vector<RobinBoundaryFace3D> robin_faces;
    // Select one material representation: spatial scalar, spatial SPD tensor,
    // or one positive isotropic coefficient per tetrahedron in mesh.elements order.
    SymmetricTensorField3D diffusivity_tensor{};
    std::vector<double> cell_diffusivity{};
    // One symmetric SPD tensor per tetrahedron in mesh-global axes.
    // Mutually exclusive with spatial scalar/tensor and cell scalar materials.
    std::vector<SymmetricTensor3D> cell_diffusivity_tensors{};
    // Optional finite volumetric source density per tetrahedron (W/m^3).
    // Additive to `source` at every quadrature point; an empty spatial callback
    // contributes zero. Negative means cooling.
    // If nonempty, follows mesh.elements order and must cover every cell.
    std::vector<double> cell_source{};
    // Stationary studies consuming absolute temperature can disallow arbitrary
    // zero-mean gauges. Every component then needs a Dirichlet or positive Robin anchor.
    bool require_anchored_components{};
};

struct ScalarDiffusion3DResult {
    // Physical field with the volume-mean convention on unanchored components.
    std::vector<double> values;
    // CG iterate and residual for its actual, temporarily pinned linear system.
    linalg::SolverResult linear_solver;
};

// Solves -div(K grad(u)) = f with continuous P1 tetrahedral elements.
// Dirichlet data may be supplied as vertex constraints or exterior faces.
// Neumann and Robin data are integrated on validated exterior triangular faces.
// Unassigned faces are insulated. Each connected component needs a Dirichlet
// value or positive Robin transfer to be anchored. Otherwise, source and
// outward-flux loads must balance; the returned solution then has zero
// volume-weighted mean on that component.
[[nodiscard]] ScalarDiffusion3DResult solve_scalar_diffusion_3d(
    const ScalarDiffusion3DProblem& problem,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::fem
