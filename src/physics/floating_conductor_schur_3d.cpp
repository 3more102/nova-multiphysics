#include "nova/physics/floating_conductor_schur_3d.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

namespace nova::physics {
namespace {

[[nodiscard]] double checked_abs_difference(const double a, const double b) {
    const long double difference =
        std::abs(static_cast<long double>(a) - static_cast<long double>(b));
    if (!std::isfinite(difference) ||
        difference > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Floating-conductor charge mismatch exceeds double range.");
    }
    return static_cast<double>(difference);
}

} // namespace

SchurFloatingConductor3DResult solve_floating_conductor_schur_3d(
    const Electrostatics3DProblem& dielectric_problem,
    const std::span<const ElectrodeBoundaryGroup3D> electrode_groups,
    const std::span<const FixedElectrodePotential3D> fixed_electrodes,
    const std::span<const FloatingElectrodeCharge3D> floating_electrodes,
    const linalg::SolverOptions& solver_options
) {
    // The Schur terminal operator describes the homogeneous (source-free)
    // linear FEM system. Impressed load vectors would introduce an affine
    // charge offset and must not silently be ignored.
    if (dielectric_problem.free_charge_density ||
        !dielectric_problem.cell_free_charge_density.empty() ||
        !dielectric_problem.displacement_flux_faces.empty()) {
        throw std::invalid_argument(
            "Floating Schur solve requires no volume free charge or "
            "prescribed normal displacement.");
    }

    SchurFloatingConductor3DResult result;
    result.discrete_capacitance = extract_discrete_capacitance_matrix_3d(
        dielectric_problem, electrode_groups, solver_options);
    result.terminals = solve_floating_electrode_potentials_3d(
        result.discrete_capacitance, fixed_electrodes, floating_electrodes);

    // Capacitance extraction has already validated unique exterior electrode
    // faces, disjoint electrode vertices, node ranges, and reference faces.
    std::set<fem::BoundaryFace3D> selected_faces;
    std::set<std::size_t> selected_nodes;
    std::vector<fem::BoundaryFace3D> all_faces;
    std::vector<int> all_ids;
    for (const auto& group : electrode_groups) {
        for (auto face : group.faces) {
            for (const std::size_t node : face) selected_nodes.insert(node);
            all_faces.push_back(face);
            all_ids.push_back(group.electrode_id);
            std::sort(face.begin(), face.end());
            selected_faces.insert(face);
        }
    }

    Electrostatics3DProblem field_problem = dielectric_problem;
    field_problem.potential_constraints.clear();
    field_problem.potential_faces.clear();

    // The basis solves ground ALL non-electrode original Dirichlet data.
    // Use identical grounding for the reconstructed spatial field; do not
    // use arbitrary original values from the input problem.
    for (const auto& constraint : dielectric_problem.potential_constraints) {
        if (!selected_nodes.contains(constraint.node)) {
            field_problem.potential_constraints.push_back({
                constraint.node, 0.0
            });
        }
    }
    const fem::ScalarField3D ground = [](double, double, double) {
        return 0.0;
    };
    for (const auto& condition : dielectric_problem.potential_faces) {
        auto face = condition.face;
        std::sort(face.begin(), face.end());
        if (!selected_faces.contains(face)) {
            field_problem.potential_faces.push_back({condition.face, ground});
        }
    }

    // Every electrode becomes a constant-potential Dirichlet conductor.
    // The floating subset was obtained by the exact P1 terminal Schur
    // complement (within iterative solve and matrix-assembly precision).
    for (const auto& group : electrode_groups) {
        const auto& ids = result.terminals.electrode_ids;
        const auto found = std::lower_bound(ids.begin(), ids.end(), group.electrode_id);
        if (found == ids.end() || *found != group.electrode_id) {
            throw std::logic_error("Schur solution has an unmatched electrode ID.");
        }
        const auto index = static_cast<std::size_t>(found - ids.begin());
        const double voltage = result.terminals.potential_volts[index];
        const auto constant_voltage = [voltage](double, double, double) {
            return voltage;
        };
        for (const auto& face : group.faces) {
            field_problem.potential_faces.push_back({face, constant_voltage});
        }
    }

    result.field = solve_electrostatics_3d(field_problem, solver_options);
    if (!result.field.linear_solver.converged) {
        throw std::runtime_error("Reconstructed floating Schur FEM field did not converge.");
    }

    // Weak-form reactions represent K*V at constrained electrode vertices.
    // They are the SAME discrete charge operator as the Schur matrix and
    // must match prescribed floating charges up to the linear solver error.
    result.discrete_reactions = evaluate_discrete_electrode_reactions_3d(
        field_problem, result.field, all_faces, all_ids);
    // Keep the legacy raw boundary estimate as a separate diagnostic.
    result.boundary_flux_charges = evaluate_electrode_charges_3d(
        field_problem, result.field, all_faces, all_ids);
    if (result.discrete_reactions.size() != result.terminals.electrode_ids.size() ||
        result.boundary_flux_charges.size() != result.discrete_reactions.size()) {
        throw std::runtime_error("Floating Schur returned incomplete electrode charges.");
    }

    // Both evaluator outputs are sorted by their unique electrode ID.
    for (std::size_t i = 0; i < result.discrete_reactions.size(); ++i) {
        const auto& reaction = result.discrete_reactions[i];
        const auto& surface = result.boundary_flux_charges[i];
        if (reaction.electrode_id != result.terminals.electrode_ids[i] ||
            surface.electrode_id != reaction.electrode_id) {
            throw std::logic_error("Floating Schur electrode charge ordering mismatch.");
        }
        result.max_flux_reaction_disagreement_coulombs = std::max(
            result.max_flux_reaction_disagreement_coulombs,
            checked_abs_difference(reaction.reaction_charge_coulombs,
                                   surface.conductor_charge_coulombs));
    }
    for (const auto& item : floating_electrodes) {
        const auto found = std::lower_bound(
            result.discrete_reactions.begin(), result.discrete_reactions.end(),
            item.electrode_id,
            [](const DiscreteElectrodeReaction3D& reaction, const int id) {
                return reaction.electrode_id < id;
            });
        if (found == result.discrete_reactions.end() ||
            found->electrode_id != item.electrode_id) {
            throw std::logic_error("Floating Schur charge constraint has no electrode.");
        }
        result.max_floating_reaction_mismatch_coulombs = std::max(
            result.max_floating_reaction_mismatch_coulombs,
            checked_abs_difference(found->reaction_charge_coulombs,
                                   item.charge_coulombs));
    }
    return result;
}

} // namespace nova::physics
