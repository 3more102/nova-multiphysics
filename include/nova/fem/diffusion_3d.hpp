#pragma once

#include "nova/fem/conductivity_3d.hpp"
#include "nova/fem/h1_space_3d.hpp"
#include "nova/linalg/conjugate_gradient.hpp"

#include <cstddef>
#include <functional>
#include <optional>
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
    ScalarField3D source;
    std::vector<DirichletConstraint3D> dirichlet_constraints;
    std::vector<DirichletBoundaryFace3D> dirichlet_faces;
    std::vector<NeumannBoundaryFace3D> neumann_faces;
    std::vector<RobinBoundaryFace3D> robin_faces;
    // Select one material representation: spatial scalar, spatial SPD tensor,
    // or one positive isotropic coefficient per tetrahedron in mesh.elements order.
    SymmetricTensorField3D diffusivity_tensor{};
    std::vector<double> cell_diffusivity{};
    // Optional finite volumetric source density per tetrahedron (W/m^3).
    // Additive to `source` at every quadrature point; negative means cooling.
    // If nonempty, follows mesh.elements order and must cover every cell.
    std::vector<double> cell_source{};
    // Opt in to a consistent weak-form global energy balance diagnostic.
    // Keeps normal solves free of additional reaction-row storage.
    bool compute_energy_balance{};
};

// Global integrated steady-state balance (SI watts for SI geometry/data).
// Outward flux is positive; generated heat is positive for a source.
// Dirichlet heat flow is calculated from the *uneliminated* FEM reactions,
// not from raw P1 element gradients at a boundary.
struct EnergyBalance3D {
    double generated_heat_rate{};
    double dirichlet_outward_heat_rate{};
    double neumann_outward_heat_rate{};
    double robin_outward_heat_rate{};
    double net_outward_heat_rate{};
    double imbalance{};           // net outgoing minus generated [W]
    double relative_imbalance{};  // normalized by max absolute physical contribution
};

struct ScalarDiffusion3DResult {
    // Physical field with the volume-mean convention on unanchored components.
    std::vector<double> values;
    // CG iterate and residual for its actual, temporarily pinned linear system.
    linalg::SolverResult linear_solver;
    // Populated only if ScalarDiffusion3DProblem::compute_energy_balance is true.
    std::optional<EnergyBalance3D> energy_balance;
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
