#include "nova/physics/electrothermal_3d.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

int main() {
    using namespace nova::fem;
    using namespace nova::physics;
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 4, 2, 2);
    Electrothermal3DProblem p;
    p.electrical.mesh = box.mesh;
    p.thermal.mesh = box.mesh;
    p.electrical.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    p.thermal.cell_diffusivity.assign(box.mesh.elements.size(), 4.0);
    for (const auto& face : box.boundary_faces) {
        const auto at_x = [&](double x) {
            return std::all_of(face.begin(), face.end(), [&](std::size_t i) {
                return box.mesh.nodes[i].x == x;
            });
        };
        if (at_x(0.0) || at_x(1.0)) {
            p.electrical.dirichlet_faces.push_back({
                face, [](double x, double, double) { return 10.0*x; }
            });
            p.thermal.dirichlet_faces.push_back({
                face, [](double, double, double) { return 300.0; }
            });
        }
    }
    const auto r = solve_one_way_electrothermal_3d(p);
    const auto range = std::minmax_element(
        r.thermal_solution.values.begin(), r.thermal_solution.values.end());
    std::cout << "NOVA 3D one-way electrothermal Joule demonstration\n"
              << "tetrahedra " << box.mesh.elements.size() << '\n'
              << "joule_power_watts " << r.total_joule_power_watts << '\n'
              << "thermal_generated_watts " << r.thermal_balance.generated_heat << '\n'
              << "heat_balance_defect_watts " << r.thermal_balance.signed_balance_defect << '\n'
              << "min_temperature_K " << *range.first << '\n'
              << "max_temperature_K " << *range.second << '\n';
    return r.electrical_solution.linear_solver.converged &&
           r.thermal_solution.linear_solver.converged &&
           std::abs(r.total_joule_power_watts - 200.0) < 1e-6 &&
           std::abs(r.thermal_balance.signed_balance_defect) < 1e-5 &&
           *range.second > 300.0 ? 0 : 1;
}
