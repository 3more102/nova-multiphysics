#include "nova/fem/diffusion_3d.hpp"
#include "nova/io/gmsh22_tetrahedral.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <utility>

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::invalid_argument(
                "Usage: nova_gmsh_diffusion_3d <mesh.msh> <dirichlet-physical-tag>"
            );
        }
        const std::string_view tag_text(argv[2]);
        int tag{};
        const auto converted = std::from_chars(
            tag_text.data(), tag_text.data() + tag_text.size(), tag
        );
        if (converted.ec != std::errc{} ||
            converted.ptr != tag_text.data() + tag_text.size() || tag <= 0) {
            throw std::invalid_argument("Dirichlet physical tag must be a positive integer.");
        }

        std::ifstream mesh_file(argv[1]);
        if (!mesh_file) {
            throw std::runtime_error("Cannot open Gmsh mesh input file.");
        }
        auto imported = nova::io::read_gmsh22_tetrahedra(mesh_file);

        nova::fem::ScalarDiffusion3DProblem problem;
        problem.mesh = std::move(imported.mesh);
        problem.diffusivity = [](double, double, double) { return 1.0; };
        problem.source = [](double, double, double) { return 1.0; };
        std::size_t anchored_faces{};
        for (std::size_t i = 0; i < imported.exterior_faces.size(); ++i) {
            if (imported.exterior_physical_tags[i] == tag) {
                problem.dirichlet_faces.push_back({
                    imported.exterior_faces[i],
                    [](double, double, double) { return 0.0; }
                });
                ++anchored_faces;
            }
        }
        if (anchored_faces == 0U) {
            throw std::invalid_argument(
                "Requested physical tag does not identify any exterior triangle."
            );
        }

        nova::linalg::SolverOptions options;
        options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
        options.relative_tolerance = 1.0e-11;
        const auto result = nova::fem::solve_scalar_diffusion_3d(problem, options);
        if (!result.linear_solver.converged) {
            throw std::runtime_error(
                "NOVA CG did not converge: " + result.linear_solver.message
            );
        }
        const auto extrema = std::minmax_element(
            result.values.begin(), result.values.end()
        );
        if (!std::isfinite(*extrema.first) || !std::isfinite(*extrema.second)) {
            throw std::runtime_error("Non-finite diffusion result.");
        }

        std::cout << "NOVA imported Gmsh P1 tetrahedral diffusion\n"
                  << "Nodes: " << problem.mesh.nodes.size() << '\n'
                  << "Tetrahedra: " << problem.mesh.elements.size() << '\n'
                  << "Dirichlet faces: " << anchored_faces << '\n'
                  << "CG iterations: " << result.linear_solver.iterations << '\n'
                  << "True residual: " << result.linear_solver.residual_norm << '\n'
                  << "Scalar range: " << *extrema.first << " to " << *extrema.second << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
