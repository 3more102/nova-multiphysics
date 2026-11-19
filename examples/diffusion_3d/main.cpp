#include "nova/fem/diffusion_3d.hpp"
#include "nova/io/vtk_legacy_3d.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string_view>

int main(int argc, char** argv) {
    using namespace nova::fem;
    const bool vtk_requested = argc == 3 &&
        std::string_view(argv[1]) == "--vtk" &&
        !std::string_view(argv[2]).empty();
    if (argc != 1 && !vtk_requested) {
        std::cerr << "Usage: nova_diffusion_3d [--vtk output.vtk]\n";
        return 1;
    }

    const auto box = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        8U, 8U, 8U
    );

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) {
        return 1.0;
    };
    problem.source = [](const double x, double, double) {
        return 1.0 + 3.0 * x;
    };

    for (const BoundaryFace3D& face : box.boundary_faces) {
        const bool on_x_min =
            problem.mesh.nodes[face[0]].x == 0.0 &&
            problem.mesh.nodes[face[1]].x == 0.0 &&
            problem.mesh.nodes[face[2]].x == 0.0;
        const bool on_x_max =
            problem.mesh.nodes[face[0]].x == 1.0 &&
            problem.mesh.nodes[face[1]].x == 1.0 &&
            problem.mesh.nodes[face[2]].x == 1.0;
        if (on_x_min || on_x_max) {
            problem.dirichlet_faces.push_back({
                face,
                [](double, double, double) { return 0.0; }
            });
        }
    }

    nova::linalg::SolverOptions solver_options;
    solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    solver_options.relative_tolerance = 1.0e-11;

    const auto result = solve_scalar_diffusion_3d(problem, solver_options);
    if (!result.linear_solver.converged) {
        std::cerr << "3D diffusion solve failed: "
                  << result.linear_solver.message << '\n';
        return 1;
    }

    double maximum_nodal_error = 0.0;
    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        const double x = problem.mesh.nodes[node].x;
        const double exact = x - 0.5 * x * x - 0.5 * x * x * x;
        maximum_nodal_error = std::max(
            maximum_nodal_error,
            std::abs(result.values[node] - exact)
        );
    }

    if (vtk_requested) {
        std::ofstream output(argv[2], std::ios::out | std::ios::trunc);
        if (!output) {
            std::cerr << "Cannot open 3D VTK output file.\n";
            return 1;
        }
        try {
            nova::io::write_legacy_vtk_scalar_3d(
                output, problem.mesh, result.values, "solution");
            output.close();
            if (!output) {
                throw std::runtime_error("Failed to finalize 3D VTK output file.");
            }
        } catch (const std::exception& error) {
            std::cerr << "3D VTK export failed: " << error.what() << '\n';
            return 1;
        }
    }

    std::cout << "NOVA 3D tetrahedral P1 diffusion\n"
              << "Nodes: " << problem.mesh.nodes.size() << '\n'
              << "Tetrahedra: " << problem.mesh.elements.size() << '\n'
              << "CG iterations: " << result.linear_solver.iterations << '\n'
              << "True residual: " << result.linear_solver.residual_norm << '\n'
              << "Maximum nodal error: " << maximum_nodal_error << '\n';
    if (vtk_requested) {
        std::cout << "VTK output: " << argv[2] << '\n';
    }
    return 0;
}
