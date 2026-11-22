#include "nova/physics/electrode_charge_3d.hpp"

#include "nova/fem/tagged_boundary_flux_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

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


std::vector<DiscreteElectrodeReaction3D> evaluate_discrete_electrode_reactions_3d(
    const Electrostatics3DProblem& problem,
    const Electrostatics3DResult& solution,
    const std::span<const fem::BoundaryFace3D> electrode_faces,
    const std::span<const int> electrode_ids
) {
    if (problem.free_charge_density ||
        !problem.cell_free_charge_density.empty() ||
        !problem.displacement_flux_faces.empty()) {
        throw std::invalid_argument(
            "Discrete electrode reactions require zero impressed volume and displacement loads.");
    }
    const auto& mesh = problem.mesh;
    mesh.validate();
    if (!solution.linear_solver.converged ||
        solution.potential.size() != mesh.nodes.size()) {
        throw std::invalid_argument(
            "Discrete electrode reactions require a converged full nodal solution.");
    }
    if (electrode_faces.empty() || electrode_faces.size() != electrode_ids.size()) {
        throw std::invalid_argument(
            "Discrete electrode reactions require nonempty equal-sized face and ID lists.");
    }
    for (double value : solution.potential) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Discrete electrode nodal voltages must be finite.");
        }
    }

    std::vector<bool> fixed(mesh.nodes.size(), false);
    for (const auto& constraint : problem.potential_constraints) {
        if (constraint.node >= fixed.size()) {
            throw std::out_of_range("Discrete reaction fixed node is outside the mesh.");
        }
        fixed[constraint.node] = true;
    }
    for (const auto& condition : problem.potential_faces) {
        for (std::size_t node : condition.face) {
            if (node >= fixed.size()) {
                throw std::out_of_range("Discrete reaction fixed face is outside the mesh.");
            }
            fixed[node] = true;
        }
    }

    // Validate actual exterior ownership independently of caller face winding.
    std::set<fem::BoundaryFace3D> exterior;
    for (auto face : fem::extract_boundary_faces_3d(mesh)) {
        std::sort(face.begin(), face.end());
        exterior.insert(face);
    }
    std::set<fem::BoundaryFace3D> visited;
    std::map<int, std::size_t> counts;
    std::map<int, double> potentials;
    std::vector<int> node_group(mesh.nodes.size(), -1);
    for (std::size_t i = 0; i < electrode_faces.size(); ++i) {
        const int id = electrode_ids[i];
        if (id < 0) {
            throw std::invalid_argument("Discrete reaction electrode IDs must be nonnegative.");
        }
        auto face = electrode_faces[i];
        for (std::size_t node : face) {
            if (node >= mesh.nodes.size()) {
                throw std::out_of_range("Discrete reaction electrode face has a missing node.");
            }
        }
        std::sort(face.begin(), face.end());
        if (face[0] == face[1] || face[1] == face[2] ||
            !exterior.contains(face) || !visited.insert(face).second) {
            throw std::invalid_argument(
                "Discrete reaction electrode faces must be distinct exterior triangles.");
        }
        ++counts[id];
        for (std::size_t node : face) {
            if (!fixed[node]) {
                throw std::invalid_argument(
                    "Discrete reaction electrode nodes must be Dirichlet constrained.");
            }
            if (node_group[node] != -1 && node_group[node] != id) {
                throw std::invalid_argument(
                    "Discrete reaction electrode groups cannot share vertices.");
            }
            node_group[node] = id;
            const double value = solution.potential[node];
            const auto [it, inserted] = potentials.try_emplace(id, value);
            if (!inserted) {
                const double scale = std::max({1.0, std::abs(value), std::abs(it->second)});
                if (std::abs(value - it->second) >
                    256.0 * std::numeric_limits<double>::epsilon() * scale) {
                    throw std::invalid_argument(
                        "Discrete reaction group nodes must be equipotential.");
                }
            }
        }
    }

    // The same material sampling as the assembled P1 volume stiffness:
    // flux_cell = -<epsilon>_cell grad(V), and K*V tested against a P1 hat
    // function is -volume * grad(N_i) dot flux_cell. No separate face normal
    // or boundary-flux extrapolation is used here.
    fem::ScalarDiffusion3DProblem weak;
    weak.mesh = mesh;
    weak.diffusivity = problem.permittivity;
    weak.diffusivity_tensor = problem.permittivity_tensor;
    weak.cell_diffusivity = problem.cell_permittivity;
    weak.cell_diffusivity_tensors = problem.cell_permittivity_tensors;
    const auto fluxes = fem::evaluate_p1_tetrahedral_flux_3d(weak, solution.potential);
    std::map<int, long double> reactions;
    for (const auto& [id, count] : counts) {
        (void)count;
        reactions.emplace(id, 0.0L);
    }
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const auto basis = fem::evaluate_tetrahedron_p1_physical_gradients_3d(
            mesh, element);
        const long double volume =
            std::abs(static_cast<long double>(basis.signed_jacobian_determinant)) / 6.0L;
        if (!(volume > 0.0L) || !std::isfinite(volume)) {
            throw std::overflow_error("Discrete electrode reaction tetra volume invalid.");
        }
        const auto& q = fluxes[element].flux;
        for (std::size_t local = 0; local < 4U; ++local) {
            const int group = node_group[mesh.elements[element][local]];
            if (group == -1) continue;
            const auto& grad = basis.gradients[local];
            const long double row_reaction = -volume * (
                static_cast<long double>(grad[0]) * q.x +
                static_cast<long double>(grad[1]) * q.y +
                static_cast<long double>(grad[2]) * q.z);
            reactions.at(group) += row_reaction;
        }
    }

    std::vector<DiscreteElectrodeReaction3D> result;
    result.reserve(counts.size());
    for (const auto& [id, face_count] : counts) {
        const long double reaction = reactions.at(id);
        if (!std::isfinite(reaction) ||
            std::abs(reaction) >
                static_cast<long double>(std::numeric_limits<double>::max())) {
            throw std::overflow_error("Discrete electrode reaction charge overflows double.");
        }
        result.push_back({id, face_count, potentials.at(id), static_cast<double>(reaction)});
    }
    return result;
}

} // namespace nova::physics
