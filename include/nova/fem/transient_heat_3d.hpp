#pragma once

#include "nova/fem/diffusion_3d.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

using ScalarField3DTime = std::function<double(double, double, double, double)>;

struct TimeDependentDirichletConstraint3D {
    std::size_t node{};
    ScalarField3DTime value;
};

struct TimeDependentDirichletBoundaryFace3D {
    BoundaryFace3D face{};
    ScalarField3DTime value;
};

struct TimeDependentNeumannBoundaryFace3D {
    BoundaryFace3D face{};
    // Positive values denote outward heat flux, matching NeumannBoundaryFace3D.
    ScalarField3DTime outward_flux;
};

struct TimeDependentRobinBoundaryFace3D {
    BoundaryFace3D face{};
    // Outward heat flux is h(t) * (T - T_ambient(t)); h must be nonnegative.
    ScalarField3DTime transfer_coefficient;
    ScalarField3DTime ambient_value;
};

// Experimental P1 tetrahedral transient heat model:
//   C(x) dT/dt - div(K(x) grad(T)) = Q(x) + Q_time(x,t).
// The optional time_dependent_source is additive to the existing static
// spatial and cell sources. Dirichlet and Neumann boundary data may also vary
// with time; material coefficients remain time-independent. Robin data may
// vary through time_dependent_robin_faces.
struct TransientHeat3DProblem {
    ScalarDiffusion3DProblem diffusion;
    // Exactly one representation: positive spatial C(x) or one positive C per
    // tetrahedron, in mesh.elements ordering. C has units J/(m^3 K).
    ScalarField3D volumetric_capacity;
    std::vector<double> cell_volumetric_capacity;
    // Exact initial nodal field; Dirichlet values are imposed at t=start_time.
    std::vector<double> initial_temperature;
    // Optional additive Q_time(x,y,z,t); BE samples t_(n+1), CN averages both endpoints.
    // Kept last among the original fields so positional aggregate initialization remains compatible.
    ScalarField3DTime time_dependent_source;
    // Optional time-varying nodal and face Dirichlet data. Face values are
    // sampled at the face's P1 vertices at each required physical time.
    std::vector<TimeDependentDirichletConstraint3D> time_dependent_dirichlet_constraints;
    std::vector<TimeDependentDirichletBoundaryFace3D> time_dependent_dirichlet_faces;
    // Optional outward Neumann flux q_n(x,y,z,t); BE samples t_(n+1), CN
    // averages both endpoints. These faces cannot overlap other face BCs.
    std::vector<TimeDependentNeumannBoundaryFace3D> time_dependent_neumann_faces;
    // Optional time-varying Robin transfer h and ambient temperature. BE uses
    // the new endpoint; CN uses both endpoints for the operator and ambient load.
    std::vector<TimeDependentRobinBoundaryFace3D> time_dependent_robin_faces;
};

enum class TransientHeat3DScheme {
    BackwardEuler,
    CrankNicolson
};

struct TransientHeat3DOptions {
    TransientHeat3DScheme scheme{TransientHeat3DScheme::BackwardEuler};
    double time_step{0.1};
    double start_time{0.0};
    std::size_t steps{1U};
    linalg::SolverOptions solver_options{};
    // Optional observer: step 0 is the Dirichlet-projected initial field;
    // steps 1..steps are each fully solved field at its physical time.
    // The vector reference remains valid only until the callback returns.
    // Exceptions propagate and abort the solve (no extra time-step solve).
    std::function<void(std::size_t, double, const std::vector<double>&)> on_step;
    // Default preserves the pre-existing per-step full fields.
    // Set false with on_step to keep CG diagnostics but not copy every field.
    bool retain_step_solutions{true};
};

// All energies are integrated over a single time step in joules.
// Neumann energy is negative for imposed outward heat; Robin energy includes
// transfer to/from ambient. Dirichlet energy is the discrete reaction at
// prescribed nodes, NOT an equilibrated boundary-face flux.
struct TransientHeatStepDiagnostics3D {
    double time{};
    double heat_content{};    // Integral C*T dV (J)
    double source_energy{};   // Volumetric source input (J)
    double neumann_energy{};  // Prescribed surface flux input (J)
    double robin_energy{};    // Net convection input (J)
    double dirichlet_energy{};// Reaction at temperature-prescribed DOFs (J)
    double balance_residual{};// Discrete heat-budget closure (J)
};

struct TransientHeat3DResult {
    std::vector<double> temperature;
    double final_time{};
    std::vector<linalg::SolverResult> linear_solvers; // solutions empty if retain_step_solutions=false
    std::vector<TransientHeatStepDiagnostics3D> diagnostics;
};

// Consistent P1 thermal mass and existing validated 3D conductivity integration.
// One time-independent system matrix is reused across all steps; CG receives a
// warm-start state. Purely insulated meshes are permitted because positive
// thermal mass removes the steady-state constant-mode singularity.
// The fixed-step API is complemented by the experimental fixed-mesh adaptive
// step-doubling API in adaptive_heat_3d.hpp; time-varying conductivity and
// heat capacity, and phase change remain unsupported.
[[nodiscard]] TransientHeat3DResult solve_transient_heat_3d(
    const TransientHeat3DProblem& problem,
    const TransientHeat3DOptions& options = {}
);

} // namespace nova::fem
