#include "nova/fem/mesh_refinement_3d.hpp"

#include <cmath>
#include <vector>

int main() {
    using namespace nova::fem;

    const auto coarse = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);

    std::vector<double> indicators(coarse.mesh.elements.size(), 0.0);
    indicators[0] = 1.0;
    const auto marked = mark_dorfler_elements_3d(indicators, 0.5);
    if (marked.size() != 1U || marked[0] != 0U) {
        return 1;
    }

    const auto refined = refine_tetrahedra_adaptively_3d(
        coarse.mesh, marked, coarse.boundary_faces);
    if (refined.mesh.elements.size() <= coarse.mesh.elements.size() ||
        refined.parent_elements.size() != refined.mesh.elements.size()) {
        return 2;
    }

    std::vector<double> coarse_values;
    coarse_values.reserve(coarse.mesh.nodes.size());
    for (const auto& node : coarse.mesh.nodes) {
        coarse_values.push_back(
            1.0 + 2.0 * node.x - 3.0 * node.y + 4.0 * node.z);
    }
    const auto fine_values = prolongate_p1_refined_3d(refined, coarse_values);
    if (fine_values.size() != refined.mesh.nodes.size()) {
        return 3;
    }
    for (std::size_t i = 0; i < fine_values.size(); ++i) {
        const auto& node = refined.mesh.nodes[i];
        const double exact = 1.0 + 2.0 * node.x - 3.0 * node.y + 4.0 * node.z;
        if (std::abs(fine_values[i] - exact) > 2.0e-13) {
            return 4;
        }
    }

    const std::vector<double> coarse_cell_values(
        coarse.mesh.elements.size(), 2.5);
    const auto fine_cell_values =
        prolongate_cellwise_constant_refined_3d(refined, coarse_cell_values);
    if (fine_cell_values.size() != refined.parent_elements.size()) {
        return 5;
    }
    for (const double value : fine_cell_values) {
        if (value != 2.5) {
            return 6;
        }
    }
    return 0;
}
