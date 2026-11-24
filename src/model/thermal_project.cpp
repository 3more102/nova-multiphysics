#include "nova/model/thermal_project.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace nova::model {
namespace {

bool valid_id(std::string_view id) {
    const auto letter = [](char c) {
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
    };
    return !id.empty() && id.size() <= 64 && letter(id.front()) &&
        std::all_of(id.begin() + 1, id.end(), [&](char c) {
            return letter(c) || (c >= '0' && c <= '9') || c == '_' || c == '-';
        });
}

void require_kelvin(double value) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument("Absolute temperature must be finite and nonnegative in kelvin.");
    }
}

void validate_boundaries(const std::map<BoxSide, ThermalBoundary>& boundaries) {
    bool anchored = false;
    for (const auto& [side, condition] : boundaries) {
        if (side < BoxSide::XMin || side > BoxSide::ZMax) {
            throw std::invalid_argument("Unknown box side.");
        }
        if (const auto* fixed = std::get_if<FixedTemperature>(&condition)) {
            require_kelvin(fixed->kelvin);
            anchored = true;
            // Opposing faces are disjoint. All other pairs share corner nodes.
            for (const auto& [other, other_condition] : boundaries) {
                const auto* previous = std::get_if<FixedTemperature>(&other_condition);
                if (other < side && static_cast<int>(other) / 2 != static_cast<int>(side) / 2 &&
                    previous && previous->kelvin != fixed->kelvin) {
                    throw std::invalid_argument("Intersecting fixed-temperature sides disagree.");
                }
            }
        } else if (const auto* flux = std::get_if<OutwardHeatFlux>(&condition)) {
            if (!std::isfinite(flux->watts_per_m2)) {
                throw std::invalid_argument("Outward heat flux must be finite.");
            }
        } else {
            const auto& film = std::get<Convection>(condition);
            if (!std::isfinite(film.film_w_m2_k) || film.film_w_m2_k <= 0.0) {
                throw std::invalid_argument("Convection coefficient must be finite and positive.");
            }
            require_kelvin(film.ambient_kelvin);
            anchored = true;
        }
    }
    if (!anchored) throw std::invalid_argument("Absolute temperature needs a fixed or convecting side.");
}

BoxSide side_of(const fem::BoundaryFace3D& face, const fem::TetrahedralMesh3D& mesh,
                const std::array<double, 3>& lengths) {
    constexpr std::array coordinates{&fem::Point3D::x, &fem::Point3D::y, &fem::Point3D::z};
    // The structured generator sets each endpoint coordinate to the exact bound.
    for (std::size_t axis = 0; axis < coordinates.size(); ++axis) {
        for (std::size_t end = 0; end < 2; ++end) {
            const double bound = end == 0 ? 0.0 : lengths[axis];
            if (std::all_of(face.begin(), face.end(), [&](std::size_t node) {
                return mesh.nodes[node].*coordinates[axis] == bound;
            })) return static_cast<BoxSide>(2 * axis + end);
        }
    }
    throw std::logic_error("Generated exterior triangle has no box-side selection.");
}

fem::ScalarField3D constant(double value) {
    return [value](double, double, double) { return value; };
}

} // namespace

void validate_thermal_project(const ThermalProject& project) {
    std::set<std::string_view> ids;
    for (const auto* id : {&project.id, &project.geometry.id, &project.mesh.id,
                           &project.material.id, &project.study_id}) {
        if (!valid_id(*id) || !ids.insert(*id).second) {
            throw std::invalid_argument("Project entity IDs must be distinct ASCII identifiers.");
        }
    }
    std::size_t tetrahedra = 6;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double length = project.geometry.lengths_m[axis];
        if (!std::isfinite(length) || length <= 0.0) {
            throw std::invalid_argument("Box lengths must be finite and positive in metres.");
        }
        const auto count = project.mesh.divisions[axis];
        if (count == 0 || count > max_project_tetrahedra / tetrahedra) {
            throw std::invalid_argument("Mesh subdivisions must be positive and within the project cell limit.");
        }
        tetrahedra *= count;
    }
    if (!std::isfinite(project.material.conductivity_w_m_k) || project.material.conductivity_w_m_k <= 0.0 ||
        !std::isfinite(project.generation_w_m3)) {
        throw std::invalid_argument("Thermal conductivity must be finite and positive; heat generation must be finite.");
    }
    const auto& options = project.solver;
    using P = linalg::PreconditionerType;
    if (options.max_iterations == 0 || options.max_iterations > max_project_iterations ||
        !std::isfinite(options.relative_tolerance) || options.relative_tolerance < 0 || options.relative_tolerance >= 1 ||
        !std::isfinite(options.absolute_tolerance) || options.absolute_tolerance < 0 ||
        (options.relative_tolerance == 0 && options.absolute_tolerance == 0) ||
        (options.preconditioner != P::None && options.preconditioner != P::Jacobi &&
         options.preconditioner != P::IncompleteCholesky0)) {
        throw std::invalid_argument("Invalid project CG settings.");
    }
    validate_boundaries(project.boundaries);
}

fem::ScalarDiffusion3DProblem prepare_thermal_project(const ThermalProject& project) {
    validate_thermal_project(project);
    const auto& l = project.geometry.lengths_m;
    const auto& n = project.mesh.divisions;
    auto box = fem::make_box_mesh_3d(0, l[0], 0, l[1], 0, l[2], n[0], n[1], n[2]);
    fem::ScalarDiffusion3DProblem problem;
    problem.mesh = std::move(box.mesh);
    problem.diffusivity = constant(project.material.conductivity_w_m_k);
    problem.source = constant(project.generation_w_m3);
    for (const auto& face : box.boundary_faces) {
        const auto it = project.boundaries.find(side_of(face, problem.mesh, l));
        if (it == project.boundaries.end()) continue;
        const auto& boundary = it->second;
        if (const auto* fixed = std::get_if<FixedTemperature>(&boundary)) {
            problem.dirichlet_faces.push_back({face, constant(fixed->kelvin)});
        } else if (const auto* flux = std::get_if<OutwardHeatFlux>(&boundary)) {
            problem.neumann_faces.push_back({face, constant(flux->watts_per_m2)});
        } else {
            const auto& film = std::get<Convection>(boundary);
            problem.robin_faces.push_back({face, constant(film.film_w_m2_k), constant(film.ambient_kelvin)});
        }
    }
    return problem;
}

ThermalStudyResult solve_thermal_project(const ThermalProject& project) {
    auto problem = prepare_thermal_project(project);
    auto temperature = fem::solve_scalar_diffusion_3d(problem, project.solver);
    if (!temperature.linear_solver.converged) {
        throw std::runtime_error("Thermal CG failed: " + temperature.linear_solver.message);
    }
    for (double kelvin : temperature.values) require_kelvin(kelvin);
    auto balance = fem::evaluate_diffusion_energy_balance_3d(problem, temperature.values);
    auto fluxes = fem::evaluate_p1_tetrahedral_flux_3d(problem, temperature.values);
    return {std::move(problem.mesh), std::move(temperature), std::move(fluxes), balance, project.generation_w_m3};
}

} // namespace nova::model
