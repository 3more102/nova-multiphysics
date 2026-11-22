#include "nova/fem/adaptive_diffusion_3d.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

[[nodiscard]] BoundaryFace3D canonical_face(BoundaryFace3D face) {
    std::sort(face.begin(), face.end());
    return face;
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

// Every exterior face omitted from the essential-boundary set has either
// the prescribed Neumann flux or the physical zero (insulation) flux.
// This prevents an unlabelled insulated face from disappearing from eta.
[[nodiscard]] DiffusionResidualIndicators3D estimate(
    const ScalarDiffusion3DProblem& problem,
    const ScalarDiffusion3DResult& solved,
    const std::vector<BoundaryFace3D>& exterior
) {
    const auto essential = dirichlet_map(problem);
    const auto neumann = neumann_map(problem);
    const ScalarField3D zero = [](double, double, double) { return 0.0; };
    const ScalarField3D source = problem.source ? problem.source : zero;
    std::vector<NeumannResidualFace3D> natural;
    natural.reserve(exterior.size());
    for (const auto& face : exterior) {
        const auto key = canonical_face(face);
        if (essential.contains(key)) continue;
        const auto found = neumann.find(key);
        natural.push_back({
            face,
            found == neumann.end() ? zero : found->second
        });
    }
    return estimate_diffusion_residual_error_3d(
        problem.mesh, solved.values, problem.cell_diffusivity, source, natural);
}

} // namespace

AdaptiveDiffusionStep3D solve_adaptive_diffusion_step_3d(
    const ScalarDiffusion3DProblem& problem,
    const double dorfler_theta,
    const linalg::SolverOptions& solver_options
) {
    // Reject unsupported inputs instead of silently solving a different PDE
    // or computing an indicator inconsistent with its material/source model.
    if (problem.diffusivity || problem.diffusivity_tensor ||
        !problem.cell_diffusivity_tensors.empty() ||
        problem.cell_diffusivity.size() != problem.mesh.elements.size() ||
        !problem.cell_source.empty() || !problem.robin_faces.empty()) {
        throw std::invalid_argument(
            "Adaptive 3D diffusion requires per-cell scalar conductivity, "
            "spatial or zero source, and Dirichlet/Neumann/insulated faces only.");
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
    const auto& refined = result.adaptation.refinement;
    result.refined_problem = problem;
    result.refined_problem.mesh = refined.mesh;
    result.refined_problem.cell_diffusivity = prolongate_cellwise_constant_refined_3d(
        refined, problem.cell_diffusivity);
    result.refined_problem.dirichlet_faces.clear();
    result.refined_problem.neumann_faces.clear();

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
        }
        // Remaining faces are naturally insulated.
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
