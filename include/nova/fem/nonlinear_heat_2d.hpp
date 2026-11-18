#pragma once

#include "nova/fem/transient_heat_2d.hpp"
#include "nova/nonlinear/newton_solver.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace nova::fem {

// Spatial temperature-dependent conductivity k(x,y,T) and its exact
// partial derivative dk/dT. NOVA deliberately does not use finite
// differences for its material tangent.
using ThermalConstitutiveField2D =
    std::function<double(double, double, double)>;

struct NonlinearHeat2DProblem {
    TriangleMesh2D mesh;
    // rho*c(x,y) > 0 and is independent of T and time.
    SpatialScalarField2D volumetric_capacity;
    ThermalConstitutiveField2D conductivity;
    ThermalConstitutiveField2D conductivity_temperature_derivative;
    TransientScalarField2D source;
    SpatialScalarField2D initial_temperature;
    // Optional exact P1/P2 DOF initial state (takes precedence over field).
    std::vector<double> initial_dof_values;
    // h(x,y,t) >= 0 on Robin boundaries; no temperature-dependent h yet.
    std::vector<TransientDirichletEdge2D> dirichlet_edges;
    std::vector<TransientNeumannEdge2D> neumann_edges;
    std::vector<TransientRobinEdge2D> robin_edges;
};

struct NonlinearHeat2DOptions {
    H1TriangleOrder order{H1TriangleOrder::P1};
    double start_time{0.0};
    double time_step{0.1};
    std::size_t steps{1U};
    nonlinear::NewtonOptions newton_options{};
};

struct NonlinearHeatSystem2D {
    // Residual of the discrete Backward Euler system and exact dR/dT.
    // Strong Dirichlet rows use R_i=T_i-g and identity Jacobian rows.
    std::vector<double> residual;
    linalg::CsrMatrix jacobian;
};

// Public assembler permits finite-difference verification of the
// exact material Jacobian at any finite trial state.
// Time is the END of the Backward Euler interval.
[[nodiscard]] NonlinearHeatSystem2D assemble_nonlinear_heat_be_system_2d(
    const NonlinearHeat2DProblem& problem,
    H1TriangleOrder order,
    std::span<const double> previous_temperature,
    std::span<const double> trial_temperature,
    double time_step,
    double time
);

struct NonlinearHeat2DResult {
    // Final P1/P2 nodal solution (vertices then P2 edge midpoints).
    std::vector<double> temperature;
    double final_time{};
    std::vector<nonlinear::NewtonResult> newton_steps;
};

// Solves rho*c*T_t - div(k(x,y,T)*grad(T)) = Q using Backward Euler,
// a consistent H1 FEM residual/Jacobian, and the shared damped Newton
// solver (GMRES or BiCGSTAB). Neumann outward flux and Robin convection
// are supported, with h independent of T. Boundary assignments cannot
// overlap on a mesh edge. Failed Newton steps throw; no silent success.
[[nodiscard]] NonlinearHeat2DResult solve_nonlinear_heat_2d(
    const NonlinearHeat2DProblem& problem,
    const NonlinearHeat2DOptions& options = {}
);

} // namespace nova::fem
