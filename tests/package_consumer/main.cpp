#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/poisson_1d.hpp"
#include "nova/fem/transient_heat_3d.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <cmath>
#include <vector>

int main() {
    // Exercise 1D FEM and CSR algebra across the installed library boundary.
    nova::fem::Poisson1DProblem p;
    p.elements = 16;
    p.source = [](double) { return 2.0; };
    const auto one_d = nova::fem::solve_poisson_1d(p);
    if (one_d.values.size() != 17 ||
        std::abs(one_d.values[8] - 0.25) > 1.0e-7) return 1;

    const auto a = nova::linalg::CsrMatrix::from_triplets(2, 2,
        {{0, 0, 2.0}, {1, 1, 3.0}});
    const std::vector<double> x{4.0, 5.0};
    const auto y = a.multiply(x);
    if (y.size() != 2 || std::abs(y[0] - 8.0) > 1.0e-12 ||
        std::abs(y[1] - 15.0) > 1.0e-12) return 2;

    // Independent installed 3D API consumer: -Delta(u)=0, u=x, with
    // fixed temperatures at x=0 and x=1 and insulated remaining faces.
    using namespace nova::fem;
    const auto box = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    for (const auto& face : box.boundary_faces) {
        const bool left = problem.mesh.nodes[face[0]].x == 0.0 &&
            problem.mesh.nodes[face[1]].x == 0.0 &&
            problem.mesh.nodes[face[2]].x == 0.0;
        const bool right = problem.mesh.nodes[face[0]].x == 1.0 &&
            problem.mesh.nodes[face[1]].x == 1.0 &&
            problem.mesh.nodes[face[2]].x == 1.0;
        if (left || right) {
            problem.dirichlet_faces.push_back({face,
                [](double px, double, double) { return px; }});
        }
    }
    const auto three_d = solve_scalar_diffusion_3d(problem);
    if (!three_d.linear_solver.converged ||
        three_d.values.size() != problem.mesh.nodes.size()) return 3;
    for (std::size_t i = 0; i < three_d.values.size(); ++i) {
        if (std::abs(three_d.values[i] - problem.mesh.nodes[i].x) > 1.0e-9)
            return 4;
    }
    // Verify the installed SDK exposes and links the new transient 3D API.
    // A uniform source on insulated, variable-material tetrahedra yields
    // dT/dt = Q/C exactly (the P1 conductivity term vanishes).
    TransientHeat3DProblem transient;
    transient.diffusion.mesh = box.mesh;
    transient.diffusion.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    transient.diffusion.source = [](double, double, double) { return 4.0; };
    transient.cell_volumetric_capacity.assign(box.mesh.elements.size(), 2.0);
    transient.initial_temperature.assign(box.mesh.nodes.size(), 300.0);
    TransientHeat3DOptions time;
    time.time_step = 0.25;
    time.steps = 2U;
    time.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    const auto heated = solve_transient_heat_3d(transient, time);
    if (heated.temperature.size() != box.mesh.nodes.size() ||
        heated.diagnostics.size() != time.steps) return 5;
    for (const double t : heated.temperature) {
        if (std::abs(t - 301.0) > 1.0e-8) return 6;
    }
    for (const auto& d : heated.diagnostics) {
        if (std::abs(d.balance_residual) > 1.0e-8) return 7;
    }

    return 0;
}
