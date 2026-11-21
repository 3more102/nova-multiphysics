#include "nova/physics/electrode_charge_3d.hpp"

#include "nova/fem/tagged_boundary_flux_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <set>
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


namespace {

using ElectrodeFaceKey3D = std::array<std::size_t, 3>;

[[nodiscard]] ElectrodeFaceKey3D canonical_electrode_face_3d(
    fem::BoundaryFace3D face
) {
    std::sort(face.begin(), face.end());
    return face;
}

} // namespace

CapacitanceMatrix3D extract_capacitance_matrix_3d(
    const Electrostatics3DProblem& dielectric_problem,
    const std::span<const ElectrodeBoundaryGroup3D> electrodes,
    const linalg::SolverOptions& solver_options
) {
    if (electrodes.empty()) {
        throw std::invalid_argument("Capacitance extraction requires at least one electrode.");
    }
    dielectric_problem.mesh.validate();

    std::vector<ElectrodeBoundaryGroup3D> groups(
        electrodes.begin(), electrodes.end()
    );
    std::sort(groups.begin(), groups.end(),
        [](const auto& a, const auto& b) {
            return a.electrode_id < b.electrode_id;
        });

    const std::size_t node_count = dielectric_problem.mesh.nodes.size();
    const std::size_t no_group = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> node_group(node_count, no_group);
    std::set<ElectrodeFaceKey3D> electrode_face_keys;
    std::vector<fem::BoundaryFace3D> all_faces;
    std::vector<int> all_ids;

    for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
        const auto& group = groups[group_index];
        if (group.electrode_id < 0) {
            throw std::invalid_argument("Electrode IDs must be nonnegative.");
        }
        if (group_index > 0U &&
            groups[group_index - 1U].electrode_id == group.electrode_id) {
            throw std::invalid_argument("Capacitance electrode IDs must be unique.");
        }
        if (group.faces.empty()) {
            throw std::invalid_argument("Every capacitance electrode needs at least one face.");
        }
        if (group.faces.size() > all_faces.max_size() - all_faces.size() ||
            group.faces.size() > all_ids.max_size() - all_ids.size()) {
            throw std::length_error("Capacitance electrode face list is too large.");
        }

        for (const auto& face : group.faces) {
            const auto key = canonical_electrode_face_3d(face);
            if (key[0] == key[1] || key[1] == key[2]) {
                throw std::invalid_argument("An electrode face must contain three distinct vertices.");
            }
            for (const std::size_t node : face) {
                if (node >= node_count) {
                    throw std::out_of_range("Capacitance electrode face references a missing vertex.");
                }
                if (node_group[node] != no_group &&
                    node_group[node] != group_index) {
                    throw std::invalid_argument(
                        "Different capacitance electrodes cannot share mesh vertices."
                    );
                }
                node_group[node] = group_index;
            }
            if (!electrode_face_keys.insert(key).second) {
                throw std::invalid_argument("Capacitance electrode faces must be unique.");
            }
            all_faces.push_back(face);
            all_ids.push_back(group.electrode_id);
        }
    }

    // Capacitance is the linear, charge-free response of the dielectric. The
    // caller's material model is retained; volume charge and prescribed
    // displacement loads are intentionally excluded.
    Electrostatics3DProblem grounded_problem = dielectric_problem;
    grounded_problem.free_charge_density = {};
    grounded_problem.cell_free_charge_density.clear();
    grounded_problem.displacement_flux_faces.clear();
    grounded_problem.potential_constraints.clear();
    grounded_problem.potential_faces.clear();

    bool has_ground_reference = false;
    for (const auto& condition : dielectric_problem.potential_constraints) {
        if (condition.node >= node_count) {
            throw std::out_of_range("Ground reference constraint references a missing vertex.");
        }
        // A constraint on an electrode vertex is replaced by the basis voltage.
        if (node_group[condition.node] != no_group) {
            continue;
        }
        grounded_problem.potential_constraints.push_back({condition.node, 0.0});
        has_ground_reference = true;
    }

    const fem::ScalarField3D zero_potential =
        [](double, double, double) { return 0.0; };
    for (const auto& condition : dielectric_problem.potential_faces) {
        for (const std::size_t node : condition.face) {
            if (node >= node_count) {
                throw std::out_of_range("Ground reference face references a missing vertex.");
            }
        }
        const auto key = canonical_electrode_face_3d(condition.face);
        if (electrode_face_keys.contains(key)) {
            // Replace an existing fixed-potential electrode face by its basis value.
            continue;
        }
        for (const std::size_t node : condition.face) {
            if (node_group[node] != no_group) {
                throw std::invalid_argument(
                    "A non-electrode grounded face cannot share vertices with an electrode."
                );
            }
        }
        grounded_problem.potential_faces.push_back({
            condition.face, zero_potential
        });
        has_ground_reference = true;
    }

    if (groups.size() == 1U && !has_ground_reference) {
        throw std::invalid_argument(
            "A single electrode requires a separate grounded potential reference."
        );
    }

    const std::size_t count = groups.size();
    CapacitanceMatrix3D result;
    result.electrode_ids.reserve(count);
    for (const auto& group : groups) {
        result.electrode_ids.push_back(group.electrode_id);
    }
    if (count > result.matrix_farads.max_size() / count) {
        throw std::length_error("Capacitance matrix dimensions overflow.");
    }
    result.matrix_farads.assign(count * count, 0.0);

    for (std::size_t basis_index = 0; basis_index < count; ++basis_index) {
        Electrostatics3DProblem basis_problem = grounded_problem;
        for (std::size_t group_index = 0; group_index < count; ++group_index) {
            const double basis_voltage = group_index == basis_index ? 1.0 : 0.0;
            const auto value = [basis_voltage](double, double, double) {
                return basis_voltage;
            };
            for (const auto& face : groups[group_index].faces) {
                basis_problem.potential_faces.push_back({face, value});
            }
        }

        const auto solved = solve_electrostatics_3d(
            basis_problem, solver_options
        );
        if (!solved.linear_solver.converged) {
            throw std::runtime_error(
                "Capacitance basis electrostatic solve did not converge."
            );
        }
        const auto charges = evaluate_electrode_charges_3d(
            basis_problem, solved, all_faces, all_ids
        );
        if (charges.size() != count) {
            throw std::runtime_error(
                "Capacitance basis solve returned an incomplete electrode charge set."
            );
        }
        for (std::size_t electrode_index = 0;
             electrode_index < count; ++electrode_index) {
            const double charge = charges[electrode_index].conductor_charge_coulombs;
            if (!std::isfinite(charge)) {
                throw std::overflow_error("Capacitance matrix contains a non-finite entry.");
            }
            result.matrix_farads[electrode_index * count + basis_index] = charge;
        }
    }

    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1U; j < count; ++j) {
            const double a = result.matrix_farads[i * count + j];
            const double b = result.matrix_farads[j * count + i];
            const double defect = std::abs(a - b);
            if (!std::isfinite(defect)) {
                throw std::overflow_error("Capacitance reciprocity defect is non-finite.");
            }
            result.max_reciprocity_defect_farads =
                std::max(result.max_reciprocity_defect_farads, defect);
            const double scale = std::max(std::abs(a), std::abs(b));
            if (scale > 0.0) {
                result.max_relative_reciprocity_defect =
                    std::max(result.max_relative_reciprocity_defect, defect / scale);
            }
        }
    }
    return result;
}

} // namespace nova::physics
