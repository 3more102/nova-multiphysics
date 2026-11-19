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

struct TransientHeat2DProblem {
    TriangleMesh2D mesh;
    // rho*c > 0; k > 0. Both are spatially variable, time independent.
    SpatialScalarField2D volumetric_capacity;
    SpatialScalarField2D conductivity;
    // Volumetric power density f(x,y,t), evaluated at t_(n+1).
    TransientScalarField2D source;
    SpatialScalarField2D initial_temperature;
    // Omitted edges have homogeneous natural (zero outward flux) conditions.
    std::vector<TransientDirichletEdge2D> dirichlet_edges;
};

struct TransientHeat2DOptions {
    H1TriangleOrder order{H1TriangleOrder::P1};
    double time_step{0.1};
    std::size_t steps{1U};
    linalg::SolverOptions solver_options{};
};

struct TransientHeat2DResult {
    // Final H1 DOF values, including P2 edge-midpoint DOFs.
    std::vector<double> temperature;
    double final_time{};
    // One linear solver diagnostic per successful time step.
    std::vector<linalg::SolverResult> linear_solvers;
};

// Implicit Backward Euler: (M + dt*K)u_(n+1) = M*u_n + dt*F(t_(n+1)).
// Prescribed boundary-edge temperatures are applied strongly and symmetrically
// at each time level. The initial field is projected onto t=0 boundary data.
// Homogeneous Neumann applies on all other edges. Natural time-dependent
// Neumann/Robin data and time-dependent material coefficients are not yet supported.
[[nodiscard]] TransientHeat2DResult solve_transient_heat_2d(
    const TransientHeat2DProblem& problem,
    const TransientHeat2DOptions& options = {}
);

} // namespace nova::fem
