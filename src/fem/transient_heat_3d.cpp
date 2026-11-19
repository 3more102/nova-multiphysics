#include "nova/fem/transient_heat_3d.hpp"

#include "heat_time_integration.hpp"
#include "p1_integration_3d.hpp"
#include "transient_heat_3d_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace nova::fem {
namespace {

struct ThermalVolume {
    ScalarH1Operators operators;
    std::vector<double> jacobians;
};

// The caller validates the mesh once before entering element loops.
[[nodiscard]] ThermalVolume assemble_volume(
    const TetrahedralMesh3D& mesh, const SpatialScalarField3D& capacity,
    const detail::ConductivityField3D& material
) {
    material.validate();
    if (!capacity) throw std::invalid_argument("Volumetric heat capacity is required.");
    std::vector<linalg::Triplet> mass, stiffness;
    if (mesh.elements.size() > mass.max_size() / 16) {
        throw std::length_error("Tetrahedral thermal matrix exceeds the supported size.");
    }
    mass.reserve(16 * mesh.elements.size());
    stiffness.reserve(16 * mesh.elements.size());
    std::vector<double> jacobians;
    jacobians.reserve(mesh.elements.size());
    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        const auto& cell = mesh.elements[e];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, e);
        const double jacobian = std::abs(basis.signed_jacobian_determinant);
        jacobians.push_back(jacobian);
        detail::P1Matrix3D local_mass{};
        const auto local_stiffness = detail::integrate_p1_volume_3d(
            mesh, cell, basis, material,
            [&](const Point3D& x, const auto& shape, double weight) {
                const double c = capacity(x.x, x.y, x.z);
                if (!std::isfinite(c) || c <= 0.0) {
                    throw std::invalid_argument("Thermal capacity must be finite and positive.");
                }
                for (std::size_t i = 0; i < 4; ++i) {
                    for (std::size_t j = 0; j < 4; ++j) {
                        local_mass[i][j] += c * weight * (shape[i] * shape[j]);
                    }
                }
            });
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = 0; j < 4; ++j) {
                mass.push_back({cell[i], cell[j], local_mass[i][j]});
                stiffness.push_back({cell[i], cell[j], local_stiffness[i][j]});
            }
        }
    }
    const auto n = mesh.nodes.size();
    return {{linalg::CsrMatrix::from_triplets(n, n, mass),
             linalg::CsrMatrix::from_triplets(n, n, stiffness)}, std::move(jacobians)};
}

[[nodiscard]] std::vector<double> boundary_values(const TransientHeat3DProblem& problem, double time) {
    std::vector<double> values(problem.mesh.nodes.size(), 0.0);
    std::vector<bool> assigned(values.size(), false);
    for (const auto& condition : problem.dirichlet_faces) {
        for (const auto i : condition.face) {
            const auto& p = problem.mesh.nodes[i];
            const double value = condition.temperature(p.x, p.y, p.z, time);
            const double scale = std::max({1.0, std::abs(values[i]), std::abs(value)});
            if (!std::isfinite(value) || (assigned[i] &&
                std::abs(value - values[i]) > 64.0 * std::numeric_limits<double>::epsilon() * scale)) {
                throw std::invalid_argument("Non-finite or conflicting 3D Dirichlet temperature.");
            }
            if (!assigned[i]) {
                values[i] = value;
                assigned[i] = true;
            }
        }
    }
    return values;
}

[[nodiscard]] detail::HeatLoads assemble_loads(
    const TransientHeat3DProblem& problem, const std::vector<double>& jacobians, double time
) {
    const auto& mesh = problem.mesh;
    const auto n = mesh.nodes.size();
    detail::HeatLoads loads;
    loads.volume_rhs.assign(n, 0.0);
    loads.neumann_rhs.assign(n, 0.0);
    loads.robin_ambient_rhs.assign(n, 0.0);
    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        const auto& cell = mesh.elements[e];
        for (const auto& q : tetrahedron_degree_two_quadrature) {
            const auto p = detail::simplex_point(mesh, cell, q.shape);
            const double source = problem.source(p.x, p.y, p.z, time);
            if (!std::isfinite(source)) {
                throw std::invalid_argument("3D heat source must be finite.");
            }
            const double weighted = source * (jacobians[e] * q.weight);
            for (std::size_t i = 0; i < 4; ++i) {
                loads.volume_rhs[cell[i]] += weighted * q.shape[i];
            }
        }
    }
    for (const auto& condition : problem.neumann_faces) {
        const auto& face = condition.face;
        const double area_scale = detail::face_jacobian(mesh, face);
        for (const auto& q : triangle_face_degree_two_quadrature) {
            const auto p = detail::simplex_point(mesh, face, q.shape);
            const double flux = condition.outward_flux(p.x, p.y, p.z, time);
            if (!std::isfinite(flux)) {
                throw std::invalid_argument("3D outward heat flux must be finite.");
            }
            for (std::size_t i = 0; i < 3; ++i) {
                loads.neumann_rhs[face[i]] -= flux * (area_scale * q.weight) * q.shape[i];
            }
        }
    }
    std::vector<linalg::Triplet> robin;
    for (const auto& condition : problem.robin_faces) {
        const auto& face = condition.face;
        const double area_scale = detail::face_jacobian(mesh, face);
        std::array<std::array<double, 3>, 3> local_robin{};
        for (const auto& q : triangle_face_degree_two_quadrature) {
            const auto p = detail::simplex_point(mesh, face, q.shape);
            const double h = condition.transfer_coefficient(p.x, p.y, p.z, time);
            const double ambient = condition.ambient_temperature(p.x, p.y, p.z, time);
            if (!std::isfinite(h) || h < 0.0 || !std::isfinite(ambient)) {
                throw std::invalid_argument("Invalid 3D Robin transfer or ambient temperature.");
            }
            const double weighted_h = h * (area_scale * q.weight);
            for (std::size_t i = 0; i < 3; ++i) {
                loads.robin_ambient_rhs[face[i]] += weighted_h * q.shape[i] * ambient;
                for (std::size_t j = 0; j < 3; ++j) {
                    local_robin[i][j] += weighted_h * (q.shape[i] * q.shape[j]);
                }
            }
        }
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                robin.push_back({face[i], face[j], local_robin[i][j]});
            }
        }
    }
    loads.robin_matrix = linalg::CsrMatrix::from_triplets(n, n, robin);
    return loads;
}

} // namespace

namespace detail {

PreparedTransientHeat3D prepare_transient_heat_3d(
    const TransientHeat3DProblem& problem
) {
    if (!problem.source ||
        (problem.initial_dof_values.empty() && !problem.initial_temperature)) {
        throw std::invalid_argument("3D heat source and initial field/DOFs are required.");
    }

    ExteriorFaceAssignments3D exterior(problem.mesh);
    std::vector<bool> constrained(problem.mesh.nodes.size(), false);
    for (const auto& bc : problem.dirichlet_faces) {
        exterior.claim(bc.face);
        if (!bc.temperature) {
            throw std::invalid_argument("3D Dirichlet temperature is required.");
        }
        for (const auto id : bc.face) {
            constrained[id] = true;
        }
    }
    for (const auto& bc : problem.neumann_faces) {
        exterior.claim(bc.face);
        if (!bc.outward_flux) {
            throw std::invalid_argument("3D Neumann flux is required.");
        }
    }
    for (const auto& bc : problem.robin_faces) {
        exterior.claim(bc.face);
        if (!bc.transfer_coefficient || !bc.ambient_temperature) {
            throw std::invalid_argument("3D Robin coefficient and ambient temperature are required.");
        }
    }

    auto volume = assemble_volume(problem.mesh, problem.volumetric_capacity,
                                  {problem.conductivity, problem.conductivity_tensor});
    return {
        std::move(volume.operators),
        std::move(volume.jacobians),
        std::move(constrained)
    };
}

std::vector<double> initial_transient_heat_state_3d(
    const TransientHeat3DProblem& problem
) {
    auto state = problem.initial_dof_values.empty()
        ? interpolate_p1_vertices_3d(problem.mesh.nodes, problem.initial_temperature)
        : problem.initial_dof_values;
    if (state.size() != problem.mesh.nodes.size() ||
        !std::all_of(state.begin(), state.end(), [](double value) { return std::isfinite(value); })) {
        throw std::invalid_argument("3D initial temperature must be finite at every P1 vertex.");
    }
    return state;
}

TransientHeat3DResult solve_prepared_transient_heat_3d(
    const TransientHeat3DProblem& problem,
    const TransientHeat3DOptions& options,
    const PreparedTransientHeat3D& prepared,
    std::vector<double> initial_state
) {
    validate_heat_time_grid(options);
    return integrate_heat(
        prepared.operators, options, std::move(initial_state), prepared.constrained,
        [&](double time) { return boundary_values(problem, time); },
        [&](double time) { return assemble_loads(problem, prepared.jacobians, time); }
    );
}

} // namespace detail

ScalarH1Operators3D assemble_scalar_h1_operators_3d(
    const TetrahedralMesh3D& mesh, const SpatialScalarField3D& capacity,
    const SpatialScalarField3D& conductivity
) {
    mesh.validate();
    return assemble_volume(mesh, capacity, {conductivity, {}}).operators;
}

ScalarH1Operators3D assemble_anisotropic_h1_operators_3d(
    const TetrahedralMesh3D& mesh, const SpatialScalarField3D& capacity,
    const SymmetricTensorField3D& conductivity_tensor
) {
    mesh.validate();
    return assemble_volume(mesh, capacity, {{}, conductivity_tensor}).operators;
}

TransientHeat3DResult solve_transient_heat_3d(
    const TransientHeat3DProblem& problem, const TransientHeat3DOptions& options
) {
    detail::validate_heat_time_grid(options);
    const auto prepared = detail::prepare_transient_heat_3d(problem);
    auto state = detail::initial_transient_heat_state_3d(problem);
    return detail::solve_prepared_transient_heat_3d(
        problem, options, prepared, std::move(state)
    );
}

} // namespace nova::fem
