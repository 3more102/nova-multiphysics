#include "nova/physics/electrostatics_3d.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace nova;
    const auto box = fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2, 2, 2
    );

    physics::Electrostatics3DProblem p;
    p.mesh = box.mesh;
    p.permittivity = [](double, double, double) { return 2.0; };
    for (const auto& face : box.boundary_faces) {
        const auto on_x = [&](double x) {
            return std::all_of(face.begin(), face.end(),
                [&](std::size_t vertex) { return p.mesh.nodes[vertex].x == x; });
        };
        if (on_x(0.0)) {
            p.potential_faces.push_back({
                face, [](double, double, double) { return 0.0; }
            });
        } else if (on_x(1.0)) {
            p.potential_faces.push_back({
                face, [](double, double, double) { return 10.0; }
            });
        }
    }

    const auto result = physics::solve_electrostatics_3d(p);
    if (!result.linear_solver.converged) {
        throw std::runtime_error("Electrostatic P1 CG solver did not converge.");
    }
    std::cout << std::setprecision(12)
              << "NOVA 3D electrostatic parallel-plate demonstration\n"
              << "Tetrahedra: " << p.mesh.elements.size() << '\n'
              << "Stored electric energy [J]: " << result.stored_energy_joules << '\n'
              << "Exact capacitor energy [J]: 100\n";
}
