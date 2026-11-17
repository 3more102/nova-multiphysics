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

struct ScalarDiffusion2DProblem {
    TriangleMesh2D mesh;
    std::function<double(double, double)> diffusivity;
    std::function<double(double, double)> source;
    std::vector<DirichletConstraint2D> dirichlet_constraints;
};

struct ScalarDiffusion2DResult {
    std::vector<double> values;
    linalg::SolverResult linear_solver;
};

// Solves -div(k grad(u)) = f on a triangular 2D mesh.
// Uses linear P1 triangles, degree-two triangle quadrature, strong
// Dirichlet elimination, and the NOVA conjugate-gradient solver.
// At least one Dirichlet constraint is required; pure Neumann problems
// and anisotropic diffusion tensors are not supported by this prototype.
[[nodiscard]] ScalarDiffusion2DResult solve_scalar_diffusion_2d(
    const ScalarDiffusion2DProblem& problem,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::fem
