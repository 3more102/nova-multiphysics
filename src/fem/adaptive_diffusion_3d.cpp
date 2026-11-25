#include "nova/fem/adaptive_diffusion_3d.hpp"

#include "nova/fem/p1_flux_3d.hpp"

#include <algorithm>
#include <cmath>
#include <array>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

[[nodiscard]] BoundaryFace3D canonical_face(BoundaryFace3D face) {
    std::sort(face.begin(), face.end());
    return face;
}

using BoundaryEdge3D = std::array<std::size_t, 2>;

[[nodiscard]] BoundaryEdge3D canonical_edge(BoundaryEdge3D edge) {
    std::sort(edge.begin(), edge.end());
    return edge;
}

[[nodiscard]] std::map<std::size_t, double> fixed_vertex_values(
    const ScalarDiffusion3DProblem& problem
) {
    std::map<std::size_t, double> values;
    for (const auto& condition : problem.dirichlet_constraints) {
        values.emplace(condition.node, condition.value);
    }
    return values;
}

[[nodiscard]] bool entirely_fixed_face(
    const BoundaryFace3D& face,
    const std::map<std::size_t, double>& values
) {
    return values.contains(face[0]) && values.contains(face[1]) &&
           values.contains(face[2]);
}

[[nodiscard]] std::map<BoundaryFace3D, ScalarField3D> dirichlet_map(
    const ScalarDiffusion3DProblem& problem
) {
    std::map<BoundaryFace3D, ScalarField3D> result;
    for (const auto& condition : problem.dirichlet_faces) {
        if (!result.emplace(canonical_face(condition.face), condition.value).second) {
            throw std::invalid_argument("Duplicate adaptive 3D Dirichlet face.");
        }
    }
    return result;
}

[[nodiscard]] std::map<BoundaryFace3D, ScalarField3D> neumann_map(
    const ScalarDiffusion3DProblem& problem
) {
    std::map<BoundaryFace3D, ScalarField3D> result;
    for (const auto& condition : problem.neumann_faces) {
        if (!result.emplace(canonical_face(condition.face), condition.outward_flux).second) {
            throw std::invalid_argument("Duplicate adaptive 3D Neumann face.");
        }
    }
    return result;
}

[[nodiscard]] std::map<BoundaryFace3D, RobinBoundaryFace3D> robin_map(
    const ScalarDiffusion3DProblem& problem
) {
    std::map<BoundaryFace3D, RobinBoundaryFace3D> result;
    for (const auto& condition : problem.robin_faces) {
        auto copy = condition;
        const auto key = canonical_face(condition.face);
        copy.face = key;
        if (!result.emplace(key, std::move(copy)).second) {
            throw std::invalid_argument("Duplicate adaptive 3D Robin face.");
        }
    }
    return result;
}

// Every exterior face omitted from the essential-boundary set has either
// prescribed Neumann/Robin data or the physical zero (insulation) flux.
[[nodiscard]] DiffusionResidualIndicators3D estimate(
    const ScalarDiffusion3DProblem& problem,
    const ScalarDiffusion3DResult& solved,
    const std::vector<BoundaryFace3D>& exterior
) {
    const auto essential = dirichlet_map(problem);
    const auto neumann = neumann_map(problem);
    const auto robin = robin_map(problem);
    const auto fixed_vertices = fixed_vertex_values(problem);
    const ScalarField3D zero = [](double, double, double) { return 0.0; };
    const ScalarField3D source = problem.source ? problem.source : zero;
    std::vector<NeumannResidualFace3D> natural;
    std::vector<RobinResidualFace3D> robin_residuals;
    natural.reserve(exterior.size());
    robin_residuals.reserve(robin.size());

    for (const auto& face : exterior) {
        const auto key = canonical_face(face);
        // A P1 face with all three vertices explicitly fixed already has
        // a prescribed trace unless a natural boundary model is explicitly
        // assigned to that face.
        if (essential.contains(key) ||
            (!neumann.contains(key) && !robin.contains(key) &&
             entirely_fixed_face(key, fixed_vertices))) {
            continue;
        }
        if (const auto neumann_condition = neumann.find(key);
            neumann_condition != neumann.end()) {
            natural.push_back({face, neumann_condition->second});
        } else if (const auto robin_condition = robin.find(key);
                   robin_condition != robin.end()) {
            robin_residuals.push_back({
                face,
                robin_condition->second.transfer_coefficient,
                robin_condition->second.ambient_value
            });
        } else {
            natural.push_back({face, zero});
        }
    }

    const auto cell_fields = evaluate_p1_tetrahedral_flux_3d(problem, solved.values);
    return estimate_diffusion_residual_error_3d(
        problem.mesh, solved.values, cell_fields, source, natural, robin_residuals,
        problem.cell_source);
}

template<class T>
[[nodiscard]] std::vector<T> transfer_cellwise_values(
    const RefinedTetrahedralMesh3D& refined,
    const std::vector<T>& coarse_values
) {
    std::vector<T> child_values;
    child_values.reserve(refined.parent_elements.size());
    for (const auto parent : refined.parent_elements) {
        if (parent >= coarse_values.size()) {
            throw std::logic_error("Adaptive 3D cell-material parent index is invalid.");
        }
        child_values.push_back(coarse_values[parent]);
    }
    return child_values;
}

[[nodiscard]] ScalarDiffusion3DProblem transfer_problem(
    const ScalarDiffusion3DProblem& problem,
    const RefinedTetrahedralMesh3D& refined,
    const std::vector<BoundaryFace3D>& exterior
) {
    const auto essential = dirichlet_map(problem);
    const auto neumann = neumann_map(problem);
    const auto robin = robin_map(problem);
    const auto fixed_vertices = fixed_vertex_values(problem);
    auto child_problem = problem;
    child_problem.mesh = refined.mesh;

    if (!problem.cell_diffusivity.empty()) {
        child_problem.cell_diffusivity =
            prolongate_cellwise_constant_refined_3d(refined, problem.cell_diffusivity);
    }
    if (!problem.cell_diffusivity_tensors.empty()) {
        child_problem.cell_diffusivity_tensors =
            transfer_cellwise_values(refined, problem.cell_diffusivity_tensors);
    }
    if (!problem.cell_source.empty()) {
        child_problem.cell_source =
            transfer_cellwise_values(refined, problem.cell_source);
    }

    child_problem.dirichlet_faces.clear();
    child_problem.neumann_faces.clear();
    child_problem.robin_faces.clear();

    for (std::size_t i = 0; i < refined.boundary_faces.size(); ++i) {
        const auto parent_index = refined.parent_boundary_faces[i];
        if (parent_index >= exterior.size()) {
            throw std::logic_error("Adaptive 3D boundary transfer parent is invalid.");
        }
        const auto& child = refined.boundary_faces[i];
        const auto parent_key = canonical_face(exterior[parent_index]);
        if (const auto dirichlet_condition = essential.find(parent_key);
            dirichlet_condition != essential.end()) {
            child_problem.dirichlet_faces.push_back({
                child, dirichlet_condition->second});
        } else if (const auto neumann_condition = neumann.find(parent_key);
                   neumann_condition != neumann.end()) {
            child_problem.neumann_faces.push_back({
                child, neumann_condition->second});
        } else if (const auto robin_condition = robin.find(parent_key);
                   robin_condition != robin.end()) {
            auto condition = robin_condition->second;
            condition.face = child;
            child_problem.robin_faces.push_back(std::move(condition));
        }
        // Remaining faces are naturally insulated.
    }

    // A boundary face fixed only through three explicit P1 vertex values
    // has a linear prescribed trace. Preserve it at new edge-midpoint DOFs;
    // otherwise those nodes would become unintended natural-boundary DOFs.
    std::set<BoundaryEdge3D> implicitly_fixed_edges;
    for (const auto& face : exterior) {
        const auto key = canonical_face(face);
        if (essential.contains(key) || neumann.contains(key) || robin.contains(key) ||
            !entirely_fixed_face(key, fixed_vertices)) {
            continue;
        }
        implicitly_fixed_edges.insert(canonical_edge({face[0], face[1]}));
        implicitly_fixed_edges.insert(canonical_edge({face[0], face[2]}));
        implicitly_fixed_edges.insert(canonical_edge({face[1], face[2]}));
    }
    for (std::size_t i = 0; i < refined.midpoint_parent_edges.size(); ++i) {
        const auto edge = canonical_edge(refined.midpoint_parent_edges[i]);
        if (!implicitly_fixed_edges.contains(edge)) continue;
        const double value = 0.5 * fixed_vertices.at(edge[0]) +
                             0.5 * fixed_vertices.at(edge[1]);
        child_problem.dirichlet_constraints.push_back({
            problem.mesh.nodes.size() + i, value});
    }

    return child_problem;
}

[[nodiscard]] ScalarDiffusion3DResult solve_checked(
    const ScalarDiffusion3DProblem& problem,
    const linalg::SolverOptions& options,
    const char* stage
) {
    auto solution = solve_scalar_diffusion_3d(problem, options);
    if (!solution.linear_solver.converged) {
        const auto& linear = solution.linear_solver;
        std::ostringstream message;
        message << stage << " adaptive diffusion 3D solve failed after "
                << linear.iterations << " iterations (residual " << linear.residual_norm
                << "). " << linear.message;
        throw std::runtime_error(message.str());
    }
    return solution;
}

} // namespace

AdaptiveDiffusionStep3D solve_adaptive_diffusion_step_3d(
    const ScalarDiffusion3DProblem& problem,
    const double dorfler_theta,
    const linalg::SolverOptions& solver_options
) {
    AdaptiveDiffusionStep3D result;
    result.coarse_solution = solve_checked(problem, solver_options, "Coarse");
    const auto exterior = extract_boundary_faces_3d(problem.mesh);
    result.coarse_indicators = estimate(problem, result.coarse_solution, exterior);
    result.adaptation = refine_tetrahedra_by_indicators_3d(
        problem.mesh, result.coarse_indicators.element_errors, dorfler_theta, exterior);
    const auto& refined = result.adaptation.refinement;
    result.transferred_p1_values = prolongate_p1_refined_3d(
        refined, result.coarse_solution.values);
    if (result.adaptation.marked_elements.empty()) {
        result.refined_problem = problem;
        result.refined_solution = result.coarse_solution;
        result.refined_indicators = result.coarse_indicators;
        return result;
    }
    result.refined_problem = transfer_problem(problem, refined, exterior);
    result.refined_solution = solve_checked(result.refined_problem, solver_options, "Refined");
    result.refined_indicators = estimate(
        result.refined_problem, result.refined_solution, refined.boundary_faces);
    return result;
}

AdaptiveDiffusionRun3D solve_adaptive_diffusion_3d(
    const ScalarDiffusion3DProblem& problem,
    const AdaptiveDiffusionOptions3D& options,
    const linalg::SolverOptions& solver_options
) {
    if (!std::isfinite(options.dorfler_theta) ||
        !(options.dorfler_theta > 0.0) || options.dorfler_theta > 1.0) {
        throw std::invalid_argument(
            "Adaptive 3D Dörfler theta must be finite and in (0, 1].");
    }
    if (!std::isfinite(options.absolute_tolerance) ||
        options.absolute_tolerance < 0.0 ||
        !std::isfinite(options.relative_tolerance) ||
        options.relative_tolerance < 0.0) {
        throw std::invalid_argument(
            "Adaptive 3D indicator tolerances must be finite and nonnegative.");
    }
    if (options.max_cycles == 0U) {
        throw std::invalid_argument(
            "Adaptive 3D multi-cycle solve requires at least one cycle.");
    }

    AdaptiveDiffusionRun3D result;
    result.final_problem = problem;
    result.final_solution = solve_checked(problem, solver_options, "Coarse");
    auto exterior = extract_boundary_faces_3d(problem.mesh);
    result.final_indicators = estimate(problem, result.final_solution, exterior);
    result.initial_indicator = result.final_indicators.total_error;
    result.target_indicator = std::max(options.absolute_tolerance,
                                       options.relative_tolerance * result.initial_indicator);
    if (!std::isfinite(result.target_indicator)) {
        throw std::overflow_error("Adaptive 3D indicator target exceeds the finite range.");
    }

    for (std::size_t cycle = 0; cycle < options.max_cycles; ++cycle) {
        const auto& mesh = result.final_problem.mesh;
        AdaptiveDiffusionCycleSummary3D summary;
        summary.cycle = cycle + 1U;
        summary.coarse_nodes = mesh.nodes.size();
        summary.coarse_elements = mesh.elements.size();
        summary.coarse_indicator = result.final_indicators.total_error;

        if (summary.coarse_indicator > result.target_indicator) {
            auto adaptation = refine_tetrahedra_by_indicators_3d(
                mesh, result.final_indicators.element_errors, options.dorfler_theta, exterior);
            summary.marked_elements = adaptation.marked_elements.size();
            if (summary.marked_elements > 0U) {
                result.final_problem = transfer_problem(
                    result.final_problem, adaptation.refinement, exterior);
                exterior = std::move(adaptation.refinement.boundary_faces);
            }
        }
        if (summary.marked_elements > 0U) {
            result.final_solution = solve_checked(result.final_problem, solver_options, "Refined");
            result.final_indicators = estimate(result.final_problem, result.final_solution, exterior);
        }

        summary.refined_nodes = result.final_problem.mesh.nodes.size();
        summary.refined_elements = result.final_problem.mesh.elements.size();
        summary.refined_indicator = result.final_indicators.total_error;
        result.cycles.push_back(summary);
        if (summary.refined_indicator <= result.target_indicator) {
            result.tolerance_satisfied = true;
            result.stop_reason = AdaptiveDiffusionStopReason3D::ToleranceSatisfied;
            break;
        }
        if (summary.marked_elements == 0U) {
            result.stop_reason = AdaptiveDiffusionStopReason3D::NoRefinementRequired;
            break;
        }
    }
    return result;
}

} // namespace nova::fem
