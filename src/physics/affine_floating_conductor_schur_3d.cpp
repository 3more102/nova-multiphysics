#include "nova/physics/affine_floating_conductor_schur_3d.hpp"

#include "nova/fem/simplex_quadrature.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace nova::physics {
namespace {

template <std::size_t N>
[[nodiscard]] fem::Point3D simplex_point(
    const fem::TetrahedralMesh3D& mesh,
    const std::array<std::size_t, N>& nodes,
    const std::array<double, N>& shape
) {
    // Identical physical-coordinate evaluation to the P1 diffusion assembler.
    const auto& origin = mesh.nodes[nodes[0]];
    fem::Point3D offset{};
    for (std::size_t i = 1; i < N; ++i) {
        const auto& point = mesh.nodes[nodes[i]];
        offset.x += (point.x - origin.x) * shape[i];
        offset.y += (point.y - origin.y) * shape[i];
        offset.z += (point.z - origin.z) * shape[i];
    }
    return {origin.x + offset.x, origin.y + offset.y, origin.z + offset.z};
}

[[nodiscard]] double face_jacobian(
    const fem::TetrahedralMesh3D& mesh, const fem::BoundaryFace3D& face
) {
    const auto& a = mesh.nodes[face[0]];
    const auto& b = mesh.nodes[face[1]];
    const auto& c = mesh.nodes[face[2]];
    const double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
    const double vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
    const double jacobian = std::hypot(
        uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
    if (!std::isfinite(jacobian) || !(jacobian > 0.0)) {
        throw std::invalid_argument("Affine Schur requires finite positive boundary face area.");
    }
    return jacobian;
}

[[nodiscard]] fem::BoundaryFace3D canonical_face(fem::BoundaryFace3D face) {
    std::sort(face.begin(), face.end());
    return face;
}

struct SelectedElectrodes {
    std::set<std::size_t> nodes;
    std::set<fem::BoundaryFace3D> faces;
    std::vector<fem::BoundaryFace3D> all_faces;
    std::vector<int> all_ids;
};

[[nodiscard]] SelectedElectrodes select_electrodes(
    std::span<const ElectrodeBoundaryGroup3D> groups
) {
    SelectedElectrodes selected;
    for (const auto& group : groups) {
        for (const auto& face : group.faces) {
            selected.faces.insert(canonical_face(face));
            selected.all_faces.push_back(face);
            selected.all_ids.push_back(group.electrode_id);
            for (std::size_t node : face) {
                selected.nodes.insert(node);
            }
        }
    }
    return selected;
}

[[nodiscard]] Electrostatics3DProblem make_terminal_problem(
    const Electrostatics3DProblem& input,
    std::span<const ElectrodeBoundaryGroup3D> groups,
    const SelectedElectrodes& selected,
    std::span<const int> sorted_ids,
    std::span<const double> voltages
) {
    if (sorted_ids.size() != voltages.size()) {
        throw std::invalid_argument("Affine Schur terminal voltage count is inconsistent.");
    }

    Electrostatics3DProblem problem = input;
    problem.potential_constraints.clear();
    problem.potential_faces.clear();

    // Match the homogeneous capacitance extraction's grounding convention.
    for (const auto& condition : input.potential_constraints) {
        if (!selected.nodes.contains(condition.node)) {
            problem.potential_constraints.push_back({condition.node, 0.0});
        }
    }
    const fem::ScalarField3D ground =
        [](double, double, double) { return 0.0; };
    for (const auto& condition : input.potential_faces) {
        if (!selected.faces.contains(canonical_face(condition.face))) {
            problem.potential_faces.push_back({condition.face, ground});
        }
    }

    for (const auto& group : groups) {
        const auto it = std::lower_bound(
            sorted_ids.begin(), sorted_ids.end(), group.electrode_id);
        if (it == sorted_ids.end() || *it != group.electrode_id) {
            throw std::logic_error("Affine Schur terminal electrode ID mismatch.");
        }
        const std::size_t index = static_cast<std::size_t>(it - sorted_ids.begin());
        const double voltage = voltages[index];
        const fem::ScalarField3D value = [voltage](double, double, double) {
            return voltage;
        };
        for (const auto& face : group.faces) {
            problem.potential_faces.push_back({face, value});
        }
    }
    return problem;
}

[[nodiscard]] std::vector<DiscreteElectrodeReaction3D> loaded_reactions(
    const Electrostatics3DProblem& problem,
    const Electrostatics3DResult& solved,
    const SelectedElectrodes& selected
) {
    // Reuse the independently validated K*V evaluator with identical material
    // quadrature. The supplied nodal field was obtained with loads; stripping
    // source descriptors here changes only the evaluator's source-free guard.
    // Below we subtract the original FEM load vector F on terminal vertices.
    auto homogeneous = problem;
    homogeneous.free_charge_density = {};
    homogeneous.cell_free_charge_density.clear();
    homogeneous.displacement_flux_faces.clear();
    auto reactions = evaluate_discrete_electrode_reactions_3d(
        homogeneous, solved, selected.all_faces, selected.all_ids);

    const auto& mesh = problem.mesh;
    const std::size_t none = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> terminal_of_node(mesh.nodes.size(), none);
    std::map<int, std::size_t> index_of_id;
    for (std::size_t i = 0; i < reactions.size(); ++i) {
        index_of_id.emplace(reactions[i].electrode_id, i);
    }
    for (std::size_t i = 0; i < selected.all_faces.size(); ++i) {
        const std::size_t terminal = index_of_id.at(selected.all_ids[i]);
        for (const auto node : selected.all_faces[i]) {
            if (terminal_of_node[node] != none &&
                terminal_of_node[node] != terminal) {
                throw std::logic_error("Affine Schur terminal nodes must be disjoint.");
            }
            terminal_of_node[node] = terminal;
        }
    }

    std::vector<long double> loads(reactions.size(), 0.0L);
    if (problem.free_charge_density || !problem.cell_free_charge_density.empty()) {
        for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
            const auto& cell = mesh.elements[e];
            const auto gradients =
                fem::evaluate_tetrahedron_p1_physical_gradients_3d(mesh, e);
            const double jacobian =
                std::abs(gradients.signed_jacobian_determinant);
            const double cell_source = problem.cell_free_charge_density.empty()
                ? 0.0 : problem.cell_free_charge_density[e];
            for (const auto& q : fem::tetrahedron_degree_two_quadrature) {
                const auto at = simplex_point(mesh, cell, q.shape);
                const double spatial = problem.free_charge_density
                    ? problem.free_charge_density(at.x, at.y, at.z) : 0.0;
                if (!std::isfinite(spatial) || !std::isfinite(cell_source)) {
                    throw std::invalid_argument("Affine Schur free charge must be finite.");
                }
                const double density = spatial + cell_source;
                if (!std::isfinite(density)) {
                    throw std::overflow_error("Affine Schur volume source is non-finite.");
                }
                const double weight = jacobian * q.weight;
                for (std::size_t i = 0; i < 4U; ++i) {
                    const auto terminal = terminal_of_node[cell[i]];
                    if (terminal != none) {
                        loads[terminal] += density * q.shape[i] * weight;
                    }
                }
            }
        }
    }

    for (const auto& condition : problem.displacement_flux_faces) {
        const double area_scale = face_jacobian(mesh, condition.face);
        for (const auto& q : fem::triangle_face_degree_two_quadrature) {
            const auto at = simplex_point(mesh, condition.face, q.shape);
            const double outward = condition.outward_flux(
                at.x, at.y, at.z);
            if (!std::isfinite(outward)) {
                throw std::invalid_argument(
                    "Affine Schur displacement flux must be finite.");
            }
            for (std::size_t i = 0; i < 3U; ++i) {
                const auto terminal = terminal_of_node[condition.face[i]];
                if (terminal != none) {
                    // Electrostatic weak RHS has - integral(D.n * N_i).
                    loads[terminal] -=
                        outward * q.shape[i] * area_scale * q.weight;
                }
            }
        }
    }

    for (std::size_t i = 0; i < reactions.size(); ++i) {
        const long double corrected =
            static_cast<long double>(reactions[i].reaction_charge_coulombs)
            - loads[i];
        if (!std::isfinite(corrected) ||
            std::abs(corrected) >
                static_cast<long double>(std::numeric_limits<double>::max())) {
            throw std::overflow_error("Affine Schur loaded reaction overflows double.");
        }
        reactions[i].reaction_charge_coulombs =
            static_cast<double>(corrected);
    }
    return reactions;
}

[[nodiscard]] double absolute_defect(double value, double expected) {
    const long double difference = std::abs(
        static_cast<long double>(value) - static_cast<long double>(expected));
    if (!std::isfinite(difference) ||
        difference > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Affine Schur charge residual exceeds double range.");
    }
    return static_cast<double>(difference);
}

} // namespace

AffineSchurFloatingConductor3DResult solve_affine_floating_conductor_schur_3d(
    const Electrostatics3DProblem& dielectric_problem,
    const std::span<const ElectrodeBoundaryGroup3D> electrode_groups,
    const std::span<const FixedElectrodePotential3D> fixed_electrodes,
    const std::span<const FloatingElectrodeCharge3D> floating_electrodes,
    const linalg::SolverOptions& solver_options
) {
    AffineSchurFloatingConductor3DResult result;
    result.discrete_capacitance = extract_discrete_capacitance_matrix_3d(
        dielectric_problem, electrode_groups, solver_options);
    const auto& ids = result.discrete_capacitance.electrode_ids;
    const auto selected = select_electrodes(electrode_groups);

    // Loaded zero-terminal-voltage solve provides the affine charge offset Q0.
    // This does not replace the homogeneous Schur basis solves.
    const std::vector<double> zeros(ids.size(), 0.0);
    auto zero_problem = make_terminal_problem(
        dielectric_problem, electrode_groups, selected, ids, zeros);
    const auto zero_field = solve_electrostatics_3d(
        zero_problem, solver_options);
    if (!zero_field.linear_solver.converged) {
        throw std::runtime_error("Affine Schur zero-voltage load solve did not converge.");
    }
    const auto zero_reactions = loaded_reactions(
        zero_problem, zero_field, selected);
    if (zero_reactions.size() != ids.size()) {
        throw std::logic_error("Affine Schur zero-bias charge count mismatch.");
    }
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (zero_reactions[i].electrode_id != ids[i]) {
            throw std::logic_error("Affine Schur zero-bias electrode order mismatch.");
        }
        result.impressed_charge_offset_coulombs.push_back(
            zero_reactions[i].reaction_charge_coulombs);
    }

    std::vector<FloatingElectrodeCharge3D> adjusted;
    adjusted.reserve(floating_electrodes.size());
    for (const auto& terminal : floating_electrodes) {
        const auto at = std::lower_bound(ids.begin(), ids.end(), terminal.electrode_id);
        if (at == ids.end() || *at != terminal.electrode_id) {
            throw std::invalid_argument("Affine Schur floating terminal ID is unknown.");
        }
        const std::size_t index = static_cast<std::size_t>(at - ids.begin());
        const double target = terminal.charge_coulombs -
            result.impressed_charge_offset_coulombs[index];
        if (!std::isfinite(target)) {
            throw std::overflow_error("Affine Schur adjusted terminal charge is non-finite.");
        }
        adjusted.push_back({terminal.electrode_id, target});
    }
    result.terminals = solve_floating_electrode_potentials_3d(
        result.discrete_capacitance, fixed_electrodes, adjusted);
    for (std::size_t i = 0; i < ids.size(); ++i) {
        result.terminals.estimated_charge_coulombs[i] +=
            result.impressed_charge_offset_coulombs[i];
    }

    const auto final_problem = make_terminal_problem(
        dielectric_problem, electrode_groups, selected,
        ids, result.terminals.potential_volts);
    result.field = solve_electrostatics_3d(final_problem, solver_options);
    if (!result.field.linear_solver.converged) {
        throw std::runtime_error("Affine Schur final field solve did not converge.");
    }
    result.loaded_discrete_reactions =
        loaded_reactions(final_problem, result.field, selected);
    result.boundary_flux_charges = evaluate_electrode_charges_3d(
        final_problem, result.field, selected.all_faces, selected.all_ids);
    if (result.loaded_discrete_reactions.size() != ids.size() ||
        result.boundary_flux_charges.size() != ids.size()) {
        throw std::logic_error("Affine Schur final electrode set is incomplete.");
    }
    for (std::size_t i = 0; i < ids.size(); ++i) {
        const auto& reaction = result.loaded_discrete_reactions[i];
        const auto& surface = result.boundary_flux_charges[i];
        if (reaction.electrode_id != ids[i] ||
            surface.electrode_id != ids[i]) {
            throw std::logic_error("Affine Schur final electrode order mismatch.");
        }
        result.max_flux_reaction_disagreement_coulombs =
            std::max(result.max_flux_reaction_disagreement_coulombs,
                absolute_defect(reaction.reaction_charge_coulombs,
                                surface.conductor_charge_coulombs));
    }
    for (const auto& terminal : floating_electrodes) {
        const auto at = std::lower_bound(ids.begin(), ids.end(), terminal.electrode_id);
        const std::size_t i = static_cast<std::size_t>(at - ids.begin());
        result.max_floating_reaction_mismatch_coulombs =
            std::max(result.max_floating_reaction_mismatch_coulombs,
                absolute_defect(
                    result.loaded_discrete_reactions[i].reaction_charge_coulombs,
                    terminal.charge_coulombs));
    }
    return result;
}

} // namespace nova::physics
