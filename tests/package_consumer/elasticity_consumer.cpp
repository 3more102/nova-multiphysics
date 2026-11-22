#include "nova/physics/linear_elasticity_3d.hpp"

#include <cmath>
#include <set>

int main() {
    using namespace nova::fem;
    using namespace nova::physics;
    const auto box = make_box_mesh_3d(0,1,0,1,0,1,2U,2U,2U);
    LinearElasticity3DProblem problem;
    problem.mesh = box.mesh;
    problem.cell_materials.assign(box.mesh.elements.size(), {100.0, 0.0});
    for (std::size_t node=0; node<box.mesh.nodes.size(); ++node) {
        if (box.mesh.nodes[node].x == 0.0) {
            problem.prescribed_displacements.push_back({node,{0,0,0}});
        }
    }
    for (const auto& face : box.boundary_faces) {
        if (box.mesh.nodes[face[0]].x == 1.0 &&
            box.mesh.nodes[face[1]].x == 1.0 &&
            box.mesh.nodes[face[2]].x == 1.0) {
            problem.face_tractions.push_back({face,{10.0,0.0,0.0}});
        }
    }
    const auto solved = solve_linear_elasticity_3d(problem);
    if (!solved.linear_solver.converged ||
        solved.displacement.size() != box.mesh.nodes.size() ||
        solved.cell_stress.size() != box.mesh.elements.size() ||
        std::abs(solved.strain_energy_joules - 0.5) > 1.0e-7) return 1;
    for (std::size_t node=0; node<box.mesh.nodes.size(); ++node) {
        if (std::abs(solved.displacement[node][0] -
                     0.1 * box.mesh.nodes[node].x) > 1.0e-8 ||
            std::abs(solved.displacement[node][1]) > 1.0e-8 ||
            std::abs(solved.displacement[node][2]) > 1.0e-8) return 2;
    }
    for (const auto& stress : solved.cell_stress) {
        if (std::abs(stress[0] - 10.0) > 1.0e-6) return 3;
    }
    return 0;
}
