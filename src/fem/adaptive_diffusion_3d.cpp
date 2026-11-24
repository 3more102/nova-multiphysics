#include "nova/fem/adaptive_diffusion_3d.hpp"

#include "nova/fem/p1_flux_3d.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <set>
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
        if (const auto found = neumann.find(key); found != neumann.end()) {
            natural.push_back({face, found->second});
        } else if (const auto found = robin.find(key); found != robin.end()) {
            robin_residuals.push_back({
                face,
                found->second.transfer_coefficient,
                found->second.ambient_value
            });
        } else {
            natural.push_back({face, zero});
        }
    }

    const auto cell_fields = evaluate_p1_tetrahedral_flux_3d(problem, solved.values);
    return estimate_diffusion_residual_error_3d(
        problem.mesh, solved.values, cell_fields, source, natural, robin_residuals);
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

} // namespace

AdaptiveDiffusionStep3D solve_adaptive_diffusion_step_3d(
    const ScalarDiffusion3DProblem& problem,
    const double dorfler_theta,
    const linalg::SolverOptions& solver_options
) {
    // Cellwise source arrays are not yet represented in the residual estimator.
    // All material forms and Robin boundaries supported by the base solver are
    // otherwise preserved through flux reconstruction and parent-cell transfer.
    if (!problem.cell_source.empty()) {
        throw std::invalid_argument(
            "Adaptive 3D diffusion currently requires a spatial or zero source; "
            "cellwise source vectors are not supported.");
    }

    AdaptiveDiffusionStep3D result;
    result.coarse_solution = solve_scalar_diffusion_3d(problem, solver_options);
    if (!result.coarse_solution.linear_solver.converged) {
        throw std::runtime_error("Coarse adaptive diffusion 3D solve did not converge.");
    }
    const auto exterior = extract_boundary_faces_3d(problem.mesh);
    result.coarse_indicators = estimate(problem, result.coarse_solution, exterior);
    result.adaptation = refine_tetrahedra_by_indicators_3d(
        problem.mesh, result.coarse_indicators.element_errors,
        dorfler_theta, exterior);

    result.transferred_p1_values = prolongate_p1_refined_3d(
        result.adaptation.refinement, result.coarse_solution.values);
    if (result.adaptation.marked_elements.empty()) {
        result.refined_problem = problem;
        result.refined_solution = result.coarse_solution;
        result.refined_indicators = result.coarse_indicators;
        return result;
    }

    const auto essential = dirichlet_map(problem);
    const auto neumann = neumann_map(problem);
    const auto robin = robin_map(problem);
    const auto fixed_vertices = fixed_vertex_values(problem);
    const auto& refined = result.adaptation.refinement;
    result.refined_problem = problem;
    result.refined_problem.mesh = refined.mesh;

    if (!problem.cell_diffusivity.empty()) {
        result.refined_problem.cell_diffusivity =
            prolongate_cellwise_constant_refined_3d(refined, problem.cell_diffusivity);
    }
    if (!problem.cell_diffusivity_tensors.empty()) {
        result.refined_problem.cell_diffusivity_tensors =
            transfer_cellwise_values(refined, problem.cell_diffusivity_tensors);
    }

    result.refined_problem.dirichlet_faces.clear();
    result.refined_problem.neumann_faces.clear();
    result.refined_problem.robin_faces.clear();

    for (std::size_t i = 0; i < refined.boundary_faces.size(); ++i) {
        const auto parent_index = refined.parent_boundary_faces[i];
        if (parent_index >= exterior.size()) {
            throw std::logic_error("Adaptive 3D boundary transfer parent is invalid.");
        }
        const auto& child = refined.boundary_faces[i];
        const auto parent_key = canonical_face(exterior[parent_index]);
        if (const auto found = essential.find(parent_key); found != essential.end()) {
            result.refined_problem.dirichlet_faces.push_back({child, found->second});
        } else if (const auto found = neumann.find(parent_key);
                   found != neumann.end()) {
            result.refined_problem.neumann_faces.push_back({child, found->second});
        } else if (const auto found = robin.find(parent_key); found != robin.end()) {
            auto condition = found->second;
            condition.face = child;
            result.refined_problem.robin_faces.push_back(std::move(condition));
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
        result.refined_problem.dirichlet_constraints.push_back({
            problem.mesh.nodes.size() + i, value});
    }

    result.refined_solution =
        solve_scalar_diffusion_3d(result.refined_problem, solver_options);
    if (!result.refined_solution.linear_solver.converged) {
        throw std::runtime_error("Refined adaptive diffusion 3D solve did not converge.");
    }
    result.refined_indicators = estimate(
        result.refined_problem, result.refined_solution, refined.boundary_faces);
    return result;
}

} // namespace nova::fem
