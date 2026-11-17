#pragma once

#include "nova/fem/mesh2d.hpp"
#include "nova/linalg/conjugate_gradient.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

struct DirichletConstraint2D {
    std::size_t node{};
    double value{};
};

struct NeumannBoundaryCondition2D {
    BoundaryEdge2D edge{};
    // Outward physical flux q_n = -k * grad(u) dot n.
    std::function<double(double, double)> outward_flux;
};

struct ScalarDiffusion2DProblem {
    TriangleMesh2D mesh;
    std::function<double(double, double)> diffusivity;
    std::function<double(double, double)> source;
    std::vector<DirichletConstraint2D> dirichlet_constraints;
    std::vector<NeumannBoundaryCondition2D> neumann_conditions;
};

struct ScalarDiffusion2DResult {
    std::vector<double> values;
    linalg::SolverResult linear_solver;
};

// Solves -div(k grad(u)) = f on a triangular 2D mesh.
// Uses linear P1 triangles, degree-four triangle quadrature, strong
// Dirichlet elimination, linear-edge Neumann integration, and CG.
// At least one Dirichlet constraint is required; pure Neumann systems
// and anisotropic diffusion tensors are not supported by this prototype.
[[nodiscard]] ScalarDiffusion2DResult solve_scalar_diffusion_2d(
    const ScalarDiffusion2DProblem& problem,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::fem
