#pragma once

#include "nova/fem/diffusion_3d.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

using ScalarField3DTime = std::function<double(double, double, double, double)>;
using SymmetricTensorField3DTime =
    std::function<SymmetricTensor3D(double, double, double, double)>;

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
    // Outward heat flux is h(x,y,z,t) * (T - T_ambient(x,y,z,t)); h >= 0.
    ScalarField3DTime transfer_coefficient;
    ScalarField3DTime ambient_value;
};

// Experimental P1 tetrahedral transient heat model. Static capacity uses
// C(x) dT/dt; when time_dependent_volumetric_capacity is supplied, storage is
// conservative: d(C(x,t) T)/dt. Conductivity may vary as K(x,y,z,t).

struct TransientHeat3DProblem {
    ScalarDiffusion3DProblem diffusion;
    // Static capacity: exactly one positive spatial C(x) or one positive C per
    // tetrahedron, in mesh.elements ordering. C has units J/(m^3 K). A dynamic
    // alternative is time_dependent_volumetric_capacity below.
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
    // Optional time-varying scalar conductivity or SPD tensor conductivity.
    // Select at most one; mutually exclusive with all static diffusion
    // conductivity representations (including cellwise values).
    ScalarField3DTime time_dependent_diffusivity;
    SymmetricTensorField3DTime time_dependent_diffusivity_tensor;
    // Optional time-varying positive volumetric heat capacity. Mutually
    // exclusive with volumetric_capacity and cell_volumetric_capacity.
    ScalarField3DTime time_dependent_volumetric_capacity;
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
    double heat_content{};    // Integral C(x,t)*T_h dV at this step's end (J)
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

// Consistent P1 tetrahedral transient heat integration. Time-varying
// conductivity and capacity are sampled at both material time levels; the
// capacity callback uses conservative d(C*T)/dt storage. Unchanged constrained
// CSR matrices and IC(0) factors are reused exactly. CG receives a warm-start
// state. The fixed-step API is complemented by the experimental fixed-mesh
// adaptive step-doubling API in adaptive_heat_3d.hpp; temperature-dependent
// constitutive laws, phase change and moving meshes remain unsupported.

[[nodiscard]] TransientHeat3DResult solve_transient_heat_3d(
    const TransientHeat3DProblem& problem,
    const TransientHeat3DOptions& options = {}
);

} // namespace nova::fem
