#include "nova/physics/electrode_charge_3d.hpp"

#include "nova/fem/tagged_boundary_flux_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace nova::physics {

std::vector<ElectrodeCharge3D> evaluate_electrode_charges_3d(
    const Electrostatics3DProblem& problem,
    const Electrostatics3DResult& solution,
    const std::span<const fem::BoundaryFace3D> electrode_faces,
    const std::span<const int> electrode_ids
) {
    const auto& mesh = problem.mesh;
    if (!solution.linear_solver.converged) {
        throw std::invalid_argument("Electrode charge requires a converged electrostatic solution.");
    }
    if (solution.potential.size() != mesh.nodes.size()) {
        throw std::invalid_argument("Electrode charge needs one potential per mesh vertex.");
    }
    if (solution.displacement.size() != mesh.elements.size()) {
        throw std::invalid_argument("Electrode charge needs one displacement per tetrahedron.");
    }
    for (const double potential : solution.potential) {
        if (!std::isfinite(potential)) {
            throw std::invalid_argument("Electrode potentials must be finite.");
        }
    }

    // Reuse the existing validated exterior-face flux integration. It finds
    // the unique owner of each triangle and constructs a normal oriented out
    // of the dielectric, independent of the caller's face winding.
    std::vector<fem::P1TetrahedralFlux3D> fields;
    fields.reserve(solution.displacement.size());
    for (const auto& displacement : solution.displacement) {
        fields.push_back({{}, displacement});
    }
    const auto integrated = fem::evaluate_tagged_boundary_flux_3d(
        mesh, electrode_faces, electrode_ids, fields
    );

    std::vector<bool> fixed(mesh.nodes.size(), false);
    for (const auto& condition : problem.potential_constraints) {
        if (condition.node >= fixed.size()) {
            throw std::invalid_argument("Prescribed electrode node is outside the mesh.");
        }
        fixed[condition.node] = true;
    }
    for (const auto& condition : problem.potential_faces) {
        for (const auto node : condition.face) {
            if (node >= fixed.size()) {
                throw std::invalid_argument("Prescribed electrode face is outside the mesh.");
            }
            fixed[node] = true;
        }
    }

    std::map<int, double> potentials;
    for (std::size_t i = 0; i < electrode_faces.size(); ++i) {
        const int id = electrode_ids[i];
        for (const auto node : electrode_faces[i]) {
            if (!fixed[node]) {
                throw std::invalid_argument(
                    "Electrode charge requires fixed-potential vertices on every selected face.");
            }
            const double v = solution.potential[node];
            const auto [it, inserted] = potentials.try_emplace(id, v);
            if (!inserted) {
                const double scale = std::max({1.0, std::abs(v), std::abs(it->second)});
                if (!(std::abs(v - it->second) <=
                      256.0 * std::numeric_limits<double>::epsilon() * scale)) {
                    throw std::invalid_argument(
                        "Every face in an electrode group must have one potential.");
                }
            }
        }
    }

    std::vector<ElectrodeCharge3D> result;
    result.reserve(integrated.size());
    for (const auto& group : integrated) {
        result.push_back({
            group.physical_tag,
            group.face_count,
            group.area,
            potentials.at(group.physical_tag),
            group.outward_heat_rate,
            -group.outward_heat_rate
        });
    }
    return result;
}

} // namespace nova::physics
