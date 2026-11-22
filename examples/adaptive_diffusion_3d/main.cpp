#include "nova/fem/adaptive_diffusion_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"
#include "nova/io/vtk_legacy_3d.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;

[[nodiscard]] bool on_x_plate(
    const nova::fem::TetrahedralMesh3D& mesh,
    const nova::fem::BoundaryFace3D& face, double x
) {
    return std::all_of(face.begin(), face.end(), [&](std::size_t node) {
        return mesh.nodes[node].x == x;
    });
}

// Independent physical check for exact u(x)=x(1-x), k=2, source=4.
// For each P1 cell, grad(u_h) is constant; integrate the squared
// gradient error exactly via first/second moments of the tetrahedron.
[[nodiscard]] double energy_error_squared(
    const nova::fem::ScalarDiffusion3DProblem& problem,
    const std::vector<double>& solution
) {
    const auto fluxes = nova::fem::evaluate_p1_tetrahedral_flux_3d(
        problem, solution);
    const auto quality = nova::fem::evaluate_mesh_quality_3d(problem.mesh);
    double total = 0.0;
    for (std::size_t i = 0; i < problem.mesh.elements.size(); ++i) {
        double sum_x = 0.0;
        double sum_x2 = 0.0;
        for (const auto node : problem.mesh.elements[i]) {
            const double x = problem.mesh.nodes[node].x;
            sum_x += x;
            sum_x2 += x * x;
        }
        const double mean_x = sum_x / 4.0;
        const double mean_x2 = (sum_x * sum_x + sum_x2) / 20.0;
        const auto& g = fluxes[i].gradient;
        const double a = 1.0 - g.x;
        const double squared = a * a - 4.0 * a * mean_x +
                               4.0 * mean_x2 + g.y * g.y + g.z * g.z;
        total += problem.cell_diffusivity[i] *
                 quality.elements[i].volume * squared;
    }
    return total;
}

void save_vtk(const fs::path& path,
              const nova::fem::TetrahedralMesh3D& mesh,
              const std::vector<double>& values) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot open VTK output: " + path.string());
    nova::io::write_legacy_vtk_scalar_3d(output, mesh, values, "temperature");
    output.close();
    if (!output) throw std::runtime_error("Cannot complete VTK output: " + path.string());
}

} // namespace

int main(int argc, char** argv) {
    try {
        const bool export_vtk = argc == 4 && std::string(argv[1]) == "--vtk";
        if (argc != 1 && !export_vtk) {
            throw std::invalid_argument(
                "Usage: nova_adaptive_diffusion_3d [--vtk coarse.vtk refined.vtk]");
        }
        fs::path coarse_path;
        fs::path refined_path;
        if (export_vtk) {
            coarse_path = fs::path(argv[2]);
            refined_path = fs::path(argv[3]);
            if (coarse_path.extension() != ".vtk" ||
                refined_path.extension() != ".vtk" ||
                coarse_path.lexically_normal() == refined_path.lexically_normal()) {
                throw std::invalid_argument(
                    "Provide two distinct .vtk output filenames.");
            }
        }

        using namespace nova::fem;
        const auto box = make_box_mesh_3d(
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
        ScalarDiffusion3DProblem problem;
        problem.mesh = box.mesh;
        problem.cell_diffusivity.assign(problem.mesh.elements.size(), 2.0);
        problem.source = [](double, double, double) { return 4.0; };
        for (const auto& face : box.boundary_faces) {
            if (on_x_plate(problem.mesh, face, 0.0) ||
                on_x_plate(problem.mesh, face, 1.0)) {
                problem.dirichlet_faces.push_back({
                    face, [](double, double, double) { return 0.0; }});
            }
        }

        nova::linalg::SolverOptions options;
        options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
        options.relative_tolerance = 1.0e-11;
        const auto step = solve_adaptive_diffusion_step_3d(
            problem, 0.65, options);
        if (!step.coarse_solution.linear_solver.converged ||
            !step.refined_solution.linear_solver.converged ||
            step.adaptation.marked_elements.empty()) {
            throw std::runtime_error("Manufactured adaptive wall solve failed.");
        }

        const double coarse_error2 =
            energy_error_squared(problem, step.coarse_solution.values);
        const double refined_error2 =
            energy_error_squared(step.refined_problem, step.refined_solution.values);
        if (!std::isfinite(coarse_error2) || !std::isfinite(refined_error2) ||
            !(coarse_error2 > 0.0) ||
            !(refined_error2 >= 0.0 && refined_error2 < coarse_error2)) {
            throw std::runtime_error(
                "Adaptive wall did not reduce the independent energy error.");
        }

        const auto coarse_quality = evaluate_mesh_quality_3d(problem.mesh);
        const auto refined_quality =
            evaluate_mesh_quality_3d(step.refined_problem.mesh);
        std::cout << std::setprecision(12)
                  << "NOVA adaptive 3D diffusion demonstration\n"
                  << "Exact solution: u(x)=x(1-x), k=2, f=4\n"
                  << "Coarse nodes/tetrahedra: " << problem.mesh.nodes.size()
                  << "/" << problem.mesh.elements.size() << '\n'
                  << "Marked tetrahedra: " << step.adaptation.marked_elements.size() << '\n'
                  << "Refined nodes/tetrahedra: "
                  << step.refined_problem.mesh.nodes.size() << "/"
                  << step.refined_problem.mesh.elements.size() << '\n'
                  << "Coarse/refined residual indicators: "
                  << step.coarse_indicators.total_error << " / "
                  << step.refined_indicators.total_error << '\n'
                  << "Coarse/refined energy-seminorm errors: "
                  << std::sqrt(coarse_error2) << " / "
                  << std::sqrt(refined_error2) << '\n'
                  << "Coarse/refined minimum mean ratios: "
                  << coarse_quality.minimum_mean_ratio << " / "
                  << refined_quality.minimum_mean_ratio << '\n';
        if (export_vtk) {
            save_vtk(coarse_path, problem.mesh, step.coarse_solution.values);
            save_vtk(refined_path, step.refined_problem.mesh,
                     step.refined_solution.values);
            std::cout << "ParaView VTK files: " << coarse_path.string()
                      << " / " << refined_path.string() << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "NOVA adaptive 3D diffusion error: "
                  << error.what() << '\n';
        return 1;
    }
}
