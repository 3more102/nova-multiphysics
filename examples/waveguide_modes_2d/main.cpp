#include "nova/fem/mesh2d.hpp"
#include "nova/physics/waveguide_modes_2d.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

int main() {
    constexpr double wavelength = 1.55e-6;
    constexpr double core_index = 1.50;
    constexpr double cladding_index = 1.45;
    constexpr double core_half_width = 1.20e-6;
    constexpr double domain_half_height = 6.0e-6;

    auto rectangle = nova::fem::make_rectangle_mesh_2d(
        -2.0e-6, 2.0e-6, -domain_half_height, domain_half_height, 2U, 40U);

    nova::physics::WaveguideModeProblem2D problem;
    problem.mesh = std::move(rectangle.mesh);
    problem.vacuum_wavelength_m = wavelength;
    problem.cladding_index = cladding_index;

    for (const std::size_t node : rectangle.boundary_nodes) {
        const double y = problem.mesh.nodes[node].y;
        if (std::abs(std::abs(y) - domain_half_height) < 1e-14) {
            problem.dirichlet_boundary_nodes.push_back(node);
        }
    }
    for (const auto& element : problem.mesh.elements) {
        const double center_y = (
            problem.mesh.nodes[element[0]].y +
            problem.mesh.nodes[element[1]].y +
            problem.mesh.nodes[element[2]].y) / 3.0;
        const double index = std::abs(center_y) < core_half_width
            ? core_index : cladding_index;
        problem.materials.push_back({index * index, 1.0});
    }

    std::cout << "NOVA 2D scalar waveguide eigenmode reference\n"
              << "P1 triangles=" << problem.mesh.elements.size()
              << " vertices=" << problem.mesh.nodes.size() << '\n'
              << std::setprecision(12);

    for (const auto polarization : {
             nova::physics::WaveguidePolarization2D::TE_Ez,
             nova::physics::WaveguidePolarization2D::TM_Hz}) {
        problem.polarization = polarization;
        const auto result = nova::physics::solve_waveguide_modes_2d(problem, 4U);
        const char* name = polarization == nova::physics::WaveguidePolarization2D::TE_Ez
            ? "TE_Ez" : "TM_Hz";
        if (result.modes.empty()) {
            std::cout << name << ": no bound modes found\n";
            continue;
        }
        const auto& fundamental = result.modes.front();
        std::cout << name << " fundamental: n_eff=" << fundamental.effective_index
                  << " beta=" << fundamental.propagation_constant_per_m
                  << " residual=" << fundamental.relative_residual
                  << " converged=" << (result.converged ? "yes" : "no") << '\n';
    }
    return 0;
}
