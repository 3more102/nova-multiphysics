#pragma once

#include "nova/fem/h1_space_2d.hpp"
#include "nova/linalg/conjugate_gradient.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

using SpatialScalarField2D = std::function<double(double, double)>;
using TransientScalarField2D = std::function<double(double, double, double)>;

// Galerkin operators for rho*c*du/dt - div(k*grad(u)) = f.
// Both matrices use continuous P1/P2 H1 ordering: vertices, then edge DOFs.
struct ScalarH1Operators2D {
    linalg::CsrMatrix mass;
    linalg::CsrMatrix stiffness;
};

[[nodiscard]] ScalarH1Operators2D assemble_scalar_h1_operators_2d(
    const TriangleMesh2D& mesh,
    H1TriangleOrder order,
    const SpatialScalarField2D& volumetric_capacity,
    const SpatialScalarField2D& conductivity
);

struct TransientDirichletEdge2D {
    BoundaryEdge2D edge{};
    TransientScalarField2D temperature;
};

struct TransientNeumannEdge2D {
    BoundaryEdge2D edge{};
    // Outward conductive heat flux q_n=-k grad(T).n, positive leaving domain.
    TransientScalarField2D outward_flux;
};

struct TransientRobinEdge2D {
    BoundaryEdge2D edge{};
    // Outward flux = h(x,y,t) * (T - ambient), with finite h >= 0.
    TransientScalarField2D transfer_coefficient;
    TransientScalarField2D ambient_temperature;
};

enum class TransientTimeScheme {
    BackwardEuler,
    CrankNicolson
};

struct TransientHeat2DProblem {
    TriangleMesh2D mesh;
    // rho*c > 0; k > 0. Both are spatially variable, time independent.
    SpatialScalarField2D volumetric_capacity;
    SpatialScalarField2D conductivity;
    // Volumetric power density f(x,y,t), evaluated at t_(n+1).
    TransientScalarField2D source;
    SpatialScalarField2D initial_temperature;
    // Optional exact H1 DOF state for restart/adaptive stepping.
    // If nonempty, overrides initial_temperature at the start time.
    std::vector<double> initial_dof_values;
    // Omitted edges have homogeneous natural (zero outward flux) conditions.
    std::vector<TransientDirichletEdge2D> dirichlet_edges;
    std::vector<TransientNeumannEdge2D> neumann_edges;
    std::vector<TransientRobinEdge2D> robin_edges;
};

struct TransientHeat2DOptions {
    H1TriangleOrder order{H1TriangleOrder::P1};
    TransientTimeScheme scheme{TransientTimeScheme::BackwardEuler};
    double time_step{0.1};
    // Physical start time; default zero preserves existing behavior.
    double start_time{0.0};
    std::size_t steps{1U};
    linalg::SolverOptions solver_options{};
};

// Per-step *integrated* thermal energy budget in the chosen physical units.
// Dirichlet energy is computed from the constrained finite-element residual.
// Energy balance includes discretization/iterative residual, not an exact
// measurement of physical flux on prescribed-temperature boundaries.
struct TransientHeatStepDiagnostics2D {
    double time{};
    double heat_content{};
    double source_energy{};
    double neumann_energy{};
    double robin_energy{};
    double dirichlet_energy{};
    double balance_residual{};
};

struct TransientHeat2DResult {
    // Final H1 DOF values, including P2 edge-midpoint DOFs.
    std::vector<double> temperature;
    double final_time{};
    // One linear solver diagnostic per successful time step.
    std::vector<linalg::SolverResult> linear_solvers;
    std::vector<TransientHeatStepDiagnostics2D> diagnostics;
};

// Theta time integration: theta=1 (Backward Euler), theta=1/2 (Crank-Nicolson).
// Time-dependent Neumann/Robin terms are integrated at both time levels.
// Boundary-edge conditions must be mutually exclusive; omitted edges are insulated.
// Prescribed temperatures are enforced symmetrically at each new time level.
// The initial field is projected onto prescribed boundary temperatures at start_time.
 // Nonempty initial_dof_values allows restarting without nodal reinterpolation.
// Material coefficients are positive, spatially varying, time-independent.
[[nodiscard]] TransientHeat2DResult solve_transient_heat_2d(
    const TransientHeat2DProblem& problem,
    const TransientHeat2DOptions& options = {}
);

} // namespace nova::fem
