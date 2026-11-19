#pragma once

#include "nova/fem/mesh2d.hpp"
#include "nova/fem/triangle_h1.hpp"
#include "nova/linalg/conjugate_gradient.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

struct DirichletConstraint2D {
    std::size_t node{};
    double value{};
};

struct DirichletBoundaryCondition2D {
    BoundaryEdge2D edge{};
    // Prescribed trace u=g. P2 also constrains the boundary edge midpoint.
    std::function<double(double, double)> value;
};

struct NeumannBoundaryCondition2D {
    BoundaryEdge2D edge{};
    // Outward physical flux q_n = -k * grad(u) dot n.
    std::function<double(double, double)> outward_flux;
};

struct RobinBoundaryCondition2D {
    BoundaryEdge2D edge{};
    // Robin convention: q_n = h * (u - u_infinity), with h >= 0.
    std::function<double(double, double)> transfer_coefficient;
    std::function<double(double, double)> ambient_value;
};

struct ScalarDiffusion2DProblem {
    TriangleMesh2D mesh;
    std::function<double(double, double)> diffusivity;
    std::function<double(double, double)> source;
    std::vector<DirichletConstraint2D> dirichlet_constraints;
    std::vector<DirichletBoundaryCondition2D> dirichlet_edges;
    std::vector<NeumannBoundaryCondition2D> neumann_conditions;
    std::vector<RobinBoundaryCondition2D> robin_conditions;
};

struct ScalarDiffusion2DResult {
    // FE field with zero integral mean on each unanchored component.
    std::vector<double> values;
    // Raw CG iterate and true residual for the symmetrically constrained system.
    linalg::SolverResult linear_solver;
};

// Solves -div(k grad(u)) = f with continuous P1 or P2 triangle elements.
// P2 results contain mesh-vertex DOFs followed by unique edge-midpoint DOFs.
// For prescribed P2 boundary traces, use dirichlet_edges so midpoint
// degrees of freedom are also constrained. Each unanchored component requires
// compatible volume and flux loads and returns a zero-area-mean field.
[[nodiscard]] ScalarDiffusion2DResult solve_scalar_diffusion_h1_2d(
    const ScalarDiffusion2DProblem& problem,
    H1TriangleOrder order,
    const linalg::SolverOptions& solver_options = {}
);

// Solves -div(k grad(u)) = f on a triangular 2D mesh.
// Uses linear P1 triangles, degree-four triangle quadrature, strong
// Dirichlet elimination, linear-edge Neumann/Robin integration, and CG.
// The P1 entry point uses the same component-wise compatibility and mean convention.
[[nodiscard]] ScalarDiffusion2DResult solve_scalar_diffusion_2d(
    const ScalarDiffusion2DProblem& problem,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::fem
