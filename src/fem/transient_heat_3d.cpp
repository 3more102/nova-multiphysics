#include "nova/fem/transient_heat_3d.hpp"

#include "nova/fem/simplex_quadrature.hpp"
#include "nova/linalg/conjugate_gradient.hpp"
#include "p1_integration_3d.hpp"
#include "transient_heat_3d_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::fem {
namespace detail {

using nova::linalg::CsrMatrix;
using nova::linalg::Triplet;

static double positive(double value, const char* what) {
    if (!std::isfinite(value) || !(value > 0.0)) {
        throw std::invalid_argument(std::string(what) + " must be finite and positive.");
    }
    return value;
}

static void impose(std::vector<std::optional<double>>& constraints,
            std::size_t node, double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("Transient 3D Dirichlet values must be finite.");
    }
    if (constraints[node]) {
        const double former = *constraints[node];
        const double scale = std::max({1.0, std::abs(former), std::abs(value)});
        if (std::abs(former - value) >
            64.0 * std::numeric_limits<double>::epsilon() * scale) {
            throw std::invalid_argument("Conflicting transient 3D Dirichlet values.");
        }
    } else {
        constraints[node] = value;
    }
}

static std::vector<std::optional<double>> evaluate_prescribed_values(
    const PreparedTransientHeat3D& prepared,
    double time
) {
    if (!std::isfinite(time) || prepared.mesh == nullptr ||
        prepared.time_dependent_prescribed.size() != prepared.prescribed.size()) {
        throw std::logic_error(
            "Transient 3D time-dependent Dirichlet preparation is inconsistent.");
    }
    auto values = prepared.prescribed;
    for (std::size_t node = 0; node < values.size(); ++node) {
        const auto& point = prepared.mesh->nodes[node];
        for (const auto& field : prepared.time_dependent_prescribed[node]) {
            if (!field) {
                throw std::logic_error(
                    "Transient 3D time-dependent Dirichlet callback is empty.");
            }
            impose(values, node, field(point.x, point.y, point.z, time));
        }
    }
    return values;
}

static double sum(const std::vector<double>& values) {
    long double total = 0.0L;
    for (double value : values) total += value;
    if (!std::isfinite(total) ||
        std::abs(total) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Transient 3D energy sum is outside the finite range.");
    }
    return static_cast<double>(total);
}

static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    long double total = 0.0L;
    for (std::size_t i = 0; i < a.size(); ++i) {
        total += static_cast<long double>(a[i]) * b[i];
    }
    if (!std::isfinite(total) ||
        std::abs(total) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Transient 3D heat content is outside the finite range.");
    }
    return static_cast<double>(total);
}

static void validate_time_grid(const TransientHeat3DOptions& options) {
    if (!std::isfinite(options.time_step) || !(options.time_step > 0.0) ||
        !std::isfinite(options.start_time) || options.steps == 0U) {
        throw std::invalid_argument("Transient 3D time step, start and step count are invalid.");
    }
    const double final_time = options.start_time +
        options.time_step * static_cast<double>(options.steps);
    if (!std::isfinite(final_time) || !(final_time > options.start_time)) {
        throw std::invalid_argument("Transient 3D final time must advance and be finite.");
    }
    double previous_time = options.start_time;
    for (std::size_t completed = 0; completed < options.steps; ++completed) {
        const double step_time = options.start_time +
            options.time_step * static_cast<double>(completed + 1U);
        if (!std::isfinite(step_time) || !(step_time > previous_time)) {
            throw std::invalid_argument(
                "Transient 3D time step is too small to advance representable physical time.");
        }
        previous_time = step_time;
    }
    if (options.scheme != TransientHeat3DScheme::BackwardEuler &&
        options.scheme != TransientHeat3DScheme::CrankNicolson) {
        throw std::invalid_argument("Unknown transient 3D time scheme.");
    }
}

PreparedTransientHeat3D prepare_transient_heat_3d(
    const TransientHeat3DProblem& problem,
    double start_time
) {
    if (!std::isfinite(start_time)) {
        throw std::invalid_argument("Transient 3D start time must be finite.");
    }
    const auto& diffusion = problem.diffusion;
    const auto& mesh = diffusion.mesh;
    // Validates geometry, connectivity, and exterior boundary assignments.
    ExteriorFaceAssignments3D exterior(mesh);
    const std::size_t n = mesh.nodes.size();
    const bool spatial_capacity = static_cast<bool>(problem.volumetric_capacity);
    const bool cell_capacity = !problem.cell_volumetric_capacity.empty();
    if (spatial_capacity == cell_capacity) {
        throw std::invalid_argument("Choose exactly one transient 3D capacity representation.");
    }
    if (cell_capacity && problem.cell_volumetric_capacity.size() != mesh.elements.size()) {
        throw std::invalid_argument("Transient 3D cell capacity length mismatch.");
    }
    if (diffusion.cell_source.size() != 0U &&
        diffusion.cell_source.size() != mesh.elements.size()) {
        throw std::invalid_argument("Transient 3D cell heat-source length mismatch.");
    }
    for (const double capacity : problem.cell_volumetric_capacity) {
        (void)positive(capacity, "Volumetric heat capacity");
    }
    for (const double source : diffusion.cell_source) {
        if (!std::isfinite(source)) {
            throw std::invalid_argument("Transient 3D cell source must be finite.");
        }
    }
    if (problem.initial_temperature.size() != n) {
        throw std::invalid_argument("Transient 3D initial temperature length mismatch.");
    }

    const DiffusionMaterials3D materials(diffusion);
    PreparedTransientHeat3D prepared;
    prepared.mesh = &mesh;
    prepared.time_dependent_source = problem.time_dependent_source;
    if (prepared.time_dependent_source) {
        prepared.element_jacobians.reserve(mesh.elements.size());
    }
    prepared.prescribed.resize(n);
    prepared.constrained.assign(n, false);
    prepared.time_dependent_prescribed.resize(n);
    std::set<std::size_t> explicit_nodes;
    for (const auto& boundary : diffusion.dirichlet_constraints) {
        if (boundary.node >= n) {
            throw std::out_of_range("Transient 3D Dirichlet node is missing.");
        }
        if (!explicit_nodes.insert(boundary.node).second) {
            throw std::invalid_argument("Duplicate transient 3D Dirichlet vertex constraint.");
        }
        impose(prepared.prescribed, boundary.node, boundary.value);
    }
    for (const auto& boundary : diffusion.dirichlet_faces) {
        exterior.claim(boundary.face);
        if (!boundary.value) {
            throw std::invalid_argument("Transient 3D Dirichlet face needs a value.");
        }
        for (const auto node : boundary.face) {
            const auto& point = mesh.nodes[node];
            impose(prepared.prescribed, node,
                   boundary.value(point.x, point.y, point.z));
        }
    }
    for (const auto& boundary : problem.time_dependent_dirichlet_constraints) {
        if (boundary.node >= n) {
            throw std::out_of_range(
                "Transient 3D time-dependent Dirichlet node is missing.");
        }
        if (!boundary.value) {
            throw std::invalid_argument(
                "Transient 3D time-dependent Dirichlet constraint needs a value.");
        }
        prepared.time_dependent_prescribed[boundary.node].push_back(boundary.value);
    }
    for (const auto& boundary : problem.time_dependent_dirichlet_faces) {
        exterior.claim(boundary.face);
        if (!boundary.value) {
            throw std::invalid_argument(
                "Transient 3D time-dependent Dirichlet face needs a value.");
        }
        for (const auto node : boundary.face) {
            prepared.time_dependent_prescribed[node].push_back(boundary.value);
        }
    }
    for (const auto& boundary : diffusion.neumann_faces) {
        exterior.claim(boundary.face);
        if (!boundary.outward_flux) {
            throw std::invalid_argument("Transient 3D Neumann face needs a flux.");
        }
    }
    for (const auto& boundary : diffusion.robin_faces) {
        exterior.claim(boundary.face);
        if (!boundary.transfer_coefficient || !boundary.ambient_value) {
            throw std::invalid_argument(
                "Transient 3D Robin face needs transfer and ambient fields.");
        }
    }

    for (const auto& boundary : problem.time_dependent_neumann_faces) {
        exterior.claim(boundary.face);
        if (!boundary.outward_flux) {
            throw std::invalid_argument(
                "Transient 3D time-dependent Neumann face needs a flux.");
        }
        const double area_scale = face_jacobian(mesh, boundary.face);
        prepared.time_dependent_neumann_faces.push_back({
            boundary.face, boundary.outward_flux, area_scale
        });
    }

    for (const auto& boundary : problem.time_dependent_robin_faces) {
        exterior.claim(boundary.face);
        if (!boundary.transfer_coefficient || !boundary.ambient_value) {
            throw std::invalid_argument(
                "Transient 3D time-dependent Robin face needs transfer and ambient fields.");
        }
        const double area_scale = face_jacobian(mesh, boundary.face);
        prepared.time_dependent_robin_faces.push_back({
            boundary.face, boundary.transfer_coefficient,
            boundary.ambient_value, area_scale
        });
    }

    for (std::size_t node = 0; node < n; ++node) {
        prepared.constrained[node] =
            prepared.prescribed[node].has_value() ||
            !prepared.time_dependent_prescribed[node].empty();
    }

    prepared.initial_state = problem.initial_temperature;
    for (const double value : prepared.initial_state) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "Transient 3D initial temperatures must be finite.");
        }
    }
    const auto initial_prescribed =
        evaluate_prescribed_values(prepared, start_time);
    for (std::size_t node = 0; node < n; ++node) {
        if (initial_prescribed[node]) {
            prepared.initial_state[node] = *initial_prescribed[node];
        }
    }

    std::vector<Triplet> mass_entries, stiffness_entries, robin_entries;
    if (mesh.elements.size() > mass_entries.max_size() / 16U ||
        mesh.elements.size() > stiffness_entries.max_size() / 16U ||
        diffusion.robin_faces.size() > robin_entries.max_size() / 27U) {
        throw std::length_error("Transient 3D matrices exceed storage limits.");
    }
    mass_entries.reserve(mesh.elements.size() * 16U);
    stiffness_entries.reserve(
        mesh.elements.size() * 16U + diffusion.robin_faces.size() * 27U);
    robin_entries.reserve(diffusion.robin_faces.size() * 27U);

    prepared.volume_load.assign(n, 0.0);
    prepared.neumann_load.assign(n, 0.0);
    prepared.ambient_load.assign(n, 0.0);
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const auto& cell = mesh.elements[element];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, element);
        if (prepared.time_dependent_source) {
            const double volume_jacobian = std::abs(basis.signed_jacobian_determinant);
            if (!std::isfinite(volume_jacobian) || !(volume_jacobian > 0.0)) {
                throw std::overflow_error("Transient 3D element volume is not representable.");
            }
            prepared.element_jacobians.push_back(volume_jacobian);
        }
        std::array<std::array<double, 4>, 4> local_mass{};
        std::array<double, 4> local_source{};
        const double extra_source =
            diffusion.cell_source.empty() ? 0.0 : diffusion.cell_source[element];
        const auto local_stiffness = integrate_p1_volume_3d(
            mesh, cell, basis, materials.on_cell(element),
            [&](const Point3D& point, const auto& shape, double weight) {
                const double capacity = cell_capacity
                    ? problem.cell_volumetric_capacity[element]
                    : positive(problem.volumetric_capacity(
                                   point.x, point.y, point.z),
                               "Volumetric heat capacity");
                const double spatial_source = diffusion.source
                    ? diffusion.source(point.x, point.y, point.z) : 0.0;
                if (!std::isfinite(spatial_source)) {
                    throw std::invalid_argument("Transient 3D source must be finite.");
                }
                const double source = spatial_source + extra_source;
                if (!std::isfinite(source)) {
                    throw std::overflow_error("Transient 3D combined source overflowed.");
                }
                for (std::size_t i = 0; i < 4; ++i) {
                    local_source[i] += source * shape[i] * weight;
                    for (std::size_t j = 0; j < 4; ++j) {
                        local_mass[i][j] += capacity * shape[i] * shape[j] * weight;
                    }
                }
            }
        );
        for (std::size_t i = 0; i < 4; ++i) {
            prepared.volume_load[cell[i]] += local_source[i];
            if (!std::isfinite(prepared.volume_load[cell[i]])) {
                throw std::overflow_error("Transient 3D source load overflowed.");
            }
            for (std::size_t j = 0; j < 4; ++j) {
                if (!std::isfinite(local_mass[i][j]) ||
                    !std::isfinite(local_stiffness[i][j])) {
                    throw std::overflow_error("Transient 3D element matrix overflowed.");
                }
                mass_entries.push_back({cell[i], cell[j], local_mass[i][j]});
                stiffness_entries.push_back({cell[i], cell[j], local_stiffness[i][j]});
            }
        }
    }

    for (const auto& boundary : diffusion.neumann_faces) {
        const double area_scale = face_jacobian(mesh, boundary.face);
        for (const auto& quadrature : triangle_face_degree_two_quadrature) {
            const auto point = simplex_point(mesh, boundary.face, quadrature.shape);
            const double flux = boundary.outward_flux(point.x, point.y, point.z);
            if (!std::isfinite(flux)) {
                throw std::invalid_argument("Transient 3D Neumann flux must be finite.");
            }
            for (std::size_t i = 0; i < 3; ++i) {
                prepared.neumann_load[boundary.face[i]] -=
                    flux * quadrature.shape[i] * area_scale * quadrature.weight;
                if (!std::isfinite(prepared.neumann_load[boundary.face[i]])) {
                    throw std::overflow_error("Transient 3D Neumann load overflowed.");
                }
            }
        }
    }
    for (const auto& boundary : diffusion.robin_faces) {
        const double area_scale = face_jacobian(mesh, boundary.face);
        for (const auto& quadrature : triangle_face_degree_two_quadrature) {
            const auto point = simplex_point(mesh, boundary.face, quadrature.shape);
            const double transfer =
                boundary.transfer_coefficient(point.x, point.y, point.z);
            const double ambient = boundary.ambient_value(point.x, point.y, point.z);
            if (!std::isfinite(transfer) || transfer < 0.0 ||
                !std::isfinite(ambient)) {
                throw std::invalid_argument(
                    "Transient 3D Robin data must be finite with h >= 0.");
            }
            const double weight = area_scale * quadrature.weight;
            for (std::size_t i = 0; i < 3; ++i) {
                prepared.ambient_load[boundary.face[i]] +=
                    transfer * ambient * quadrature.shape[i] * weight;
                if (!std::isfinite(prepared.ambient_load[boundary.face[i]])) {
                    throw std::overflow_error("Transient 3D Robin load overflowed.");
                }
                for (std::size_t j = 0; j < 3; ++j) {
                    const double value =
                        transfer * quadrature.shape[i] * quadrature.shape[j] * weight;
                    if (!std::isfinite(value)) {
                        throw std::overflow_error("Transient 3D Robin matrix overflowed.");
                    }
                    robin_entries.push_back(
                        {boundary.face[i], boundary.face[j], value});
                    stiffness_entries.push_back(
                        {boundary.face[i], boundary.face[j], value});
                }
            }
        }
    }

    prepared.total_load.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        prepared.total_load[i] =
            prepared.volume_load[i] + prepared.neumann_load[i] +
            prepared.ambient_load[i];
        if (!std::isfinite(prepared.total_load[i])) {
            throw std::overflow_error("Transient 3D total load overflowed.");
        }
    }

    prepared.mass = CsrMatrix::from_triplets(n, n, mass_entries);
    prepared.stiffness = CsrMatrix::from_triplets(n, n, stiffness_entries);
    prepared.robin = CsrMatrix::from_triplets(n, n, robin_entries);
    prepared.capacity_weights = prepared.mass.multiply(
        std::vector<double>(n, 1.0));
    return prepared;
}

namespace {

struct ConstrainedStepMatrix3D {
    CsrMatrix system;
    std::vector<Triplet> boundary_couplings;
};

ConstrainedStepMatrix3D constrain_step_matrix(
    const PreparedTransientHeat3D& prepared,
    const CsrMatrix& lhs
) {
    const std::size_t n = prepared.initial_state.size();
    std::vector<Triplet> free_entries;
    std::vector<Triplet> boundary_couplings;
    if (n > free_entries.max_size() ||
        lhs.non_zeros() > free_entries.max_size() - n) {
        throw std::length_error("Transient 3D constrained matrix is too large.");
    }
    free_entries.reserve(lhs.non_zeros() + n);
    lhs.for_each_nonzero(
        [&](std::size_t row, std::size_t column, double value) {
            if (!prepared.constrained[row] && !prepared.constrained[column]) {
                free_entries.push_back({row, column, value});
            } else if (!prepared.constrained[row] &&
                       prepared.constrained[column]) {
                boundary_couplings.push_back({row, column, value});
            }
        });
    for (std::size_t i = 0; i < n; ++i) {
        if (prepared.constrained[i]) {
            free_entries.push_back({i, i, 1.0});
        }
    }
    return {
        CsrMatrix::from_triplets(n, n, free_entries),
        std::move(boundary_couplings)
    };
}

void prepare_step_operators(
    const PreparedTransientHeat3D& prepared,
    double time_step,
    TransientHeat3DScheme scheme,
    TransientHeat3DStepCache& cache
) {
    if (cache.operators && cache.time_step == time_step &&
        cache.scheme == scheme) {
        return;
    }

    const double theta =
        scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 0.5;
    TransientHeat3DStepOperators next;
    next.lhs = CsrMatrix::linear_combination(
        prepared.mass, 1.0, prepared.stiffness, theta * time_step);
    next.rhs_operator = CsrMatrix::linear_combination(
        prepared.mass, 1.0, prepared.stiffness, -(1.0 - theta) * time_step);

    auto constrained = constrain_step_matrix(prepared, next.lhs);
    next.system = std::move(constrained.system);
    next.boundary_couplings = std::move(constrained.boundary_couplings);

    if (cache.operator_builds == std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error("Transient 3D operator-build counter overflowed.");
    }
    cache.operators = std::move(next);
    cache.time_step = time_step;
    cache.scheme = scheme;
    ++cache.operator_builds;
}

std::vector<double> evaluate_time_dependent_source_load(
    const PreparedTransientHeat3D& prepared,
    double time
) {
    if (!prepared.time_dependent_source) return {};
    if (!std::isfinite(time) || prepared.mesh == nullptr ||
        prepared.element_jacobians.size() != prepared.mesh->elements.size()) {
        throw std::logic_error("Transient 3D time-dependent source preparation is inconsistent.");
    }

    const auto& mesh = *prepared.mesh;
    std::vector<double> load(prepared.initial_state.size(), 0.0);
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const auto& cell = mesh.elements[element];
        const double jacobian = prepared.element_jacobians[element];
        for (const auto& quadrature : tetrahedron_degree_two_quadrature) {
            const auto point = simplex_point(mesh, cell, quadrature.shape);
            const double source = prepared.time_dependent_source(
                point.x, point.y, point.z, time);
            if (!std::isfinite(source)) {
                throw std::invalid_argument("Transient 3D time-dependent source must be finite.");
            }
            const double weight = jacobian * quadrature.weight;
            for (std::size_t i = 0; i < 4; ++i) {
                const double contribution = source * quadrature.shape[i] * weight;
                if (!std::isfinite(contribution)) {
                    throw std::overflow_error("Transient 3D time-dependent source load overflowed.");
                }
                load[cell[i]] += contribution;
                if (!std::isfinite(load[cell[i]])) {
                    throw std::overflow_error("Transient 3D time-dependent source accumulation overflowed.");
                }
            }
        }
    }
    return load;
}


std::vector<double> evaluate_time_dependent_neumann_load(
    const PreparedTransientHeat3D& prepared,
    double time
) {
    if (prepared.time_dependent_neumann_faces.empty()) return {};
    if (!std::isfinite(time) || prepared.mesh == nullptr) {
        throw std::logic_error(
            "Transient 3D time-dependent Neumann preparation is inconsistent.");
    }

    const auto& mesh = *prepared.mesh;
    std::vector<double> load(prepared.initial_state.size(), 0.0);
    for (const auto& boundary : prepared.time_dependent_neumann_faces) {
        for (const auto& quadrature : triangle_face_degree_two_quadrature) {
            const auto point = simplex_point(mesh, boundary.face, quadrature.shape);
            const double flux = boundary.outward_flux(
                point.x, point.y, point.z, time);
            if (!std::isfinite(flux)) {
                throw std::invalid_argument(
                    "Transient 3D time-dependent Neumann flux must be finite.");
            }
            for (std::size_t i = 0; i < 3; ++i) {
                const double contribution =
                    -flux * quadrature.shape[i] * boundary.area_scale *
                    quadrature.weight;
                if (!std::isfinite(contribution)) {
                    throw std::overflow_error(
                        "Transient 3D time-dependent Neumann load overflowed.");
                }
                load[boundary.face[i]] += contribution;
                if (!std::isfinite(load[boundary.face[i]])) {
                    throw std::overflow_error(
                        "Transient 3D time-dependent Neumann accumulation overflowed.");
                }
            }
        }
    }
    return load;
}

struct TimeDependentRobinData3D {
    CsrMatrix matrix;
    std::vector<double> ambient_load;
};

TimeDependentRobinData3D evaluate_time_dependent_robin(
    const PreparedTransientHeat3D& prepared,
    double time
) {
    if (!std::isfinite(time) || prepared.mesh == nullptr ||
        prepared.time_dependent_robin_faces.empty()) {
        throw std::logic_error(
            "Transient 3D time-dependent Robin preparation is inconsistent.");
    }
    const auto& mesh = *prepared.mesh;
    const std::size_t n = prepared.initial_state.size();
    std::vector<Triplet> matrix_entries;
    if (prepared.time_dependent_robin_faces.size() >
        matrix_entries.max_size() / 27U) {
        throw std::length_error("Transient 3D Robin matrix exceeds storage limits.");
    }
    matrix_entries.reserve(prepared.time_dependent_robin_faces.size() * 27U);
    std::vector<double> ambient_load(n, 0.0);

    for (const auto& boundary : prepared.time_dependent_robin_faces) {
        for (const auto& quadrature : triangle_face_degree_two_quadrature) {
            const auto point = simplex_point(mesh, boundary.face, quadrature.shape);
            const double transfer = boundary.transfer_coefficient(
                point.x, point.y, point.z, time);
            const double ambient = boundary.ambient_value(
                point.x, point.y, point.z, time);
            if (!std::isfinite(transfer) || transfer < 0.0 ||
                !std::isfinite(ambient)) {
                throw std::invalid_argument(
                    "Transient 3D time-dependent Robin data must be finite with h >= 0.");
            }
            const double weighted_ambient = transfer * ambient;
            const double weight = boundary.area_scale * quadrature.weight;
            if (!std::isfinite(weighted_ambient) || !std::isfinite(weight)) {
                throw std::overflow_error(
                    "Transient 3D time-dependent Robin quadrature overflowed.");
            }
            for (std::size_t i = 0; i < 3; ++i) {
                const double ambient_contribution =
                    weighted_ambient * quadrature.shape[i] * weight;
                if (!std::isfinite(ambient_contribution)) {
                    throw std::overflow_error(
                        "Transient 3D time-dependent Robin load overflowed.");
                }
                ambient_load[boundary.face[i]] += ambient_contribution;
                if (!std::isfinite(ambient_load[boundary.face[i]])) {
                    throw std::overflow_error(
                        "Transient 3D time-dependent Robin load accumulation overflowed.");
                }
                for (std::size_t j = 0; j < 3; ++j) {
                    const double value =
                        transfer * quadrature.shape[i] * quadrature.shape[j] * weight;
                    if (!std::isfinite(value)) {
                        throw std::overflow_error(
                            "Transient 3D time-dependent Robin matrix overflowed.");
                    }
                    matrix_entries.push_back({boundary.face[i], boundary.face[j], value});
                }
            }
        }
    }

    return {
        CsrMatrix::from_triplets(n, n, matrix_entries),
        std::move(ambient_load)
    };
}

} // namespace

TransientHeat3DResult solve_prepared_transient_heat_3d(
    const PreparedTransientHeat3D& prepared,
    const TransientHeat3DOptions& options,
    std::vector<double> state,
    TransientHeat3DStepCache& cache
) {
    validate_time_grid(options);
    const std::size_t n = prepared.initial_state.size();
    if (state.size() != n) {
        throw std::invalid_argument("Transient 3D initial temperature length mismatch.");
    }
    for (const double value : state) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Transient 3D initial temperatures must be finite.");
        }
    }
    const auto initial_prescribed =
        evaluate_prescribed_values(prepared, options.start_time);
    for (std::size_t node = 0; node < n; ++node) {
        if (initial_prescribed[node]) {
            state[node] = *initial_prescribed[node];
        }
    }

    prepare_step_operators(prepared, options.time_step, options.scheme, cache);
    const auto& operators = *cache.operators;
    const double dt = options.time_step;
    const double theta =
        options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 0.5;
    std::vector<double> previous_time_source;
    if (prepared.time_dependent_source && theta < 1.0) {
        previous_time_source = evaluate_time_dependent_source_load(
            prepared, options.start_time);
    }
    std::vector<double> previous_time_neumann;
    if (!prepared.time_dependent_neumann_faces.empty() && theta < 1.0) {
        previous_time_neumann = evaluate_time_dependent_neumann_load(
            prepared, options.start_time);
    }
    std::optional<TimeDependentRobinData3D> previous_time_robin;
    if (!prepared.time_dependent_robin_faces.empty() && theta < 1.0) {
        previous_time_robin = evaluate_time_dependent_robin(
            prepared, options.start_time);
    }

    TransientHeat3DResult result;
    result.final_time = options.start_time +
        dt * static_cast<double>(options.steps);
    result.linear_solvers.reserve(options.steps);
    result.diagnostics.reserve(options.steps);

    // Publish the constraint-projected initial field as frame zero.
    if (options.on_step) {
        options.on_step(0U, options.start_time, state);
    }

    for (std::size_t step = 1; step <= options.steps; ++step) {
        const double step_time = options.start_time + dt * static_cast<double>(step);
        auto next_time_source = prepared.time_dependent_source
            ? evaluate_time_dependent_source_load(prepared, step_time)
            : std::vector<double>{};
        auto next_time_neumann =
            !prepared.time_dependent_neumann_faces.empty()
                ? evaluate_time_dependent_neumann_load(prepared, step_time)
                : std::vector<double>{};
        std::optional<TimeDependentRobinData3D> next_time_robin;
        if (!prepared.time_dependent_robin_faces.empty()) {
            next_time_robin = evaluate_time_dependent_robin(prepared, step_time);
        }
        std::vector<double> weighted_time_neumann;
        if (!prepared.time_dependent_neumann_faces.empty()) {
            weighted_time_neumann.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                weighted_time_neumann[i] = theta * next_time_neumann[i];
                if (theta < 1.0) {
                    weighted_time_neumann[i] +=
                        (1.0 - theta) * previous_time_neumann[i];
                }
                if (!std::isfinite(weighted_time_neumann[i])) {
                    throw std::overflow_error(
                        "Transient 3D weighted Neumann load overflowed.");
                }
            }
        }
        const auto current_prescribed =
            evaluate_prescribed_values(prepared, step_time);

        // Dynamic Robin transfer changes the endpoint system matrix. For CN,
        // its previous-time matrix also contributes to the explicit RHS.
        std::optional<CsrMatrix> time_lhs;
        std::optional<CsrMatrix> time_rhs_operator;
        std::optional<CsrMatrix> time_system;
        std::vector<Triplet> time_boundary_couplings;
        const CsrMatrix* lhs_operator = &operators.lhs;
        const CsrMatrix* rhs_operator = &operators.rhs_operator;
        const CsrMatrix* system_operator = &operators.system;
        const std::vector<Triplet>* boundary_couplings =
            &operators.boundary_couplings;
        if (next_time_robin) {
            time_lhs.emplace(CsrMatrix::linear_combination(
                operators.lhs, 1.0, next_time_robin->matrix, theta * dt));
            lhs_operator = &*time_lhs;
            if (theta < 1.0) {
                if (!previous_time_robin) {
                    throw std::logic_error(
                        "Transient 3D previous-time Robin data are missing.");
                }
                time_rhs_operator.emplace(CsrMatrix::linear_combination(
                    operators.rhs_operator, 1.0,
                    previous_time_robin->matrix, -(1.0 - theta) * dt));
                rhs_operator = &*time_rhs_operator;
            }
            auto constrained = constrain_step_matrix(prepared, *lhs_operator);
            time_system.emplace(std::move(constrained.system));
            time_boundary_couplings =
                std::move(constrained.boundary_couplings);
            system_operator = &*time_system;
            boundary_couplings = &time_boundary_couplings;
        }

        std::vector<double> weighted_time_source;
        if (prepared.time_dependent_source) {
            weighted_time_source.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                weighted_time_source[i] = theta * next_time_source[i];
                if (theta < 1.0) {
                    weighted_time_source[i] += (1.0 - theta) * previous_time_source[i];
                }
                if (!std::isfinite(weighted_time_source[i])) {
                    throw std::overflow_error("Transient 3D weighted source load overflowed.");
                }
            }
        }
        std::vector<double> weighted_time_robin_ambient;
        if (next_time_robin) {
            weighted_time_robin_ambient.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                weighted_time_robin_ambient[i] =
                    theta * next_time_robin->ambient_load[i];
                if (theta < 1.0) {
                    weighted_time_robin_ambient[i] +=
                        (1.0 - theta) * previous_time_robin->ambient_load[i];
                }
                if (!std::isfinite(weighted_time_robin_ambient[i])) {
                    throw std::overflow_error(
                        "Transient 3D weighted Robin ambient load overflowed.");
                }
            }
        }
        std::vector<double> rhs = rhs_operator->multiply(state);
        for (std::size_t i = 0; i < n; ++i) {
            const double load = prepared.total_load[i] +
                (prepared.time_dependent_source ? weighted_time_source[i] : 0.0) +
                (!prepared.time_dependent_neumann_faces.empty()
                    ? weighted_time_neumann[i] : 0.0) +
                (next_time_robin ? weighted_time_robin_ambient[i] : 0.0);
            rhs[i] += dt * load;
            if (!std::isfinite(rhs[i])) {
                throw std::overflow_error("Transient 3D RHS overflowed.");
            }
        }
        // Preserve the assembled, unconstrained rhs for reaction diagnostics.
        const auto physical_rhs = rhs;
        std::vector<double> boundary_shift(n, 0.0);
        for (const auto& edge : *boundary_couplings) {
            if (!current_prescribed[edge.column]) {
                throw std::logic_error(
                    "Transient 3D constrained DOF has no prescribed value.");
            }
            boundary_shift[edge.row] -=
                edge.value * *current_prescribed[edge.column];
            if (!std::isfinite(boundary_shift[edge.row])) {
                throw std::overflow_error(
                    "Transient 3D boundary shift overflowed.");
            }
        }

        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] += boundary_shift[i];
            if (current_prescribed[i]) {
                rhs[i] = *current_prescribed[i];
            }
            if (!std::isfinite(rhs[i])) {
                throw std::overflow_error("Transient 3D boundary RHS overflowed.");
            }
        }

        // Solve for the small increment rather than the large absolute
        // temperature. Defect correction avoids loss of relative precision
        // when the absolute temperature is much larger than its increment.
        const auto current_product = system_operator->multiply(state);
        std::vector<double> increment_rhs(n);
        for (std::size_t i = 0; i < n; ++i) {
            increment_rhs[i] = rhs[i] - current_product[i];
            if (!std::isfinite(increment_rhs[i])) {
                throw std::overflow_error("Transient 3D correction load overflowed.");
            }
        }
        auto solved = cache.cached_cg.solve(
            *system_operator, increment_rhs, options.solver_options);
        if (!solved.converged) {
            throw std::runtime_error(
                "Transient 3D step " + std::to_string(step) +
                " failed: " + solved.message);
        }

        std::vector<double> next(n);
        for (std::size_t i = 0; i < n; ++i) {
            next[i] = state[i] + solved.solution[i];
            if (current_prescribed[i]) {
                next[i] = *current_prescribed[i];
            }
            if (!std::isfinite(next[i])) {
                throw std::overflow_error(
                    "Transient 3D corrected temperature overflowed.");
            }
        }
        if (options.retain_step_solutions) {
            solved.solution = next;
        } else {
            std::vector<double>{}.swap(solved.solution);
        }

        const auto unconstrained_lhs = lhs_operator->multiply(next);
        const auto previous_robin = prepared.robin.multiply(state);
        const auto next_robin = prepared.robin.multiply(next);
        const auto previous_time_robin_product =
            previous_time_robin
                ? previous_time_robin->matrix.multiply(state)
                : std::vector<double>{};
        const auto next_time_robin_product =
            next_time_robin
                ? next_time_robin->matrix.multiply(next)
                : std::vector<double>{};

        TransientHeatStepDiagnostics3D diagnostic;
        diagnostic.time = step_time;
        if (!std::isfinite(diagnostic.time)) {
            throw std::overflow_error("Transient 3D time overflowed.");
        }
        diagnostic.heat_content = dot(prepared.capacity_weights, next);
        diagnostic.source_energy = dt * (
            sum(prepared.volume_load) + sum(weighted_time_source));
        if (!std::isfinite(diagnostic.source_energy)) {
            throw std::overflow_error("Transient 3D source energy overflowed.");
        }
        diagnostic.neumann_energy = dt * (
            sum(prepared.neumann_load) + sum(weighted_time_neumann));
        const double robin_ambient_rate =
            sum(prepared.ambient_load) +
            sum(weighted_time_robin_ambient);
        const double robin_outward_rate =
            theta * (sum(next_robin) + sum(next_time_robin_product)) +
            (1.0 - theta) *
                (sum(previous_robin) + sum(previous_time_robin_product));
        diagnostic.robin_energy = dt *
            (robin_ambient_rate - robin_outward_rate);
        for (std::size_t i = 0; i < n; ++i) {
            if (current_prescribed[i]) {
                diagnostic.dirichlet_energy +=
                    unconstrained_lhs[i] - physical_rhs[i];
            }
        }
        diagnostic.balance_residual =
            diagnostic.heat_content - dot(prepared.capacity_weights, state) -
            diagnostic.source_energy - diagnostic.neumann_energy -
            diagnostic.robin_energy - diagnostic.dirichlet_energy;
        if (!std::isfinite(diagnostic.dirichlet_energy) ||
            !std::isfinite(diagnostic.balance_residual)) {
            throw std::overflow_error("Transient 3D energy diagnostic overflowed.");
        }

        if (theta < 1.0 && prepared.time_dependent_source) {
            previous_time_source = std::move(next_time_source);
        }
        if (theta < 1.0 && !prepared.time_dependent_neumann_faces.empty()) {
            previous_time_neumann = std::move(next_time_neumann);
        }
        if (theta < 1.0 && next_time_robin) {
            previous_time_robin = std::move(next_time_robin);
        }
        state = std::move(next);
        result.diagnostics.push_back(diagnostic);
        result.linear_solvers.push_back(std::move(solved));
        if (options.on_step) {
            options.on_step(step, diagnostic.time, state);
        }
    }

    result.temperature = std::move(state);
    return result;
}

} // namespace detail

TransientHeat3DResult solve_transient_heat_3d(
    const TransientHeat3DProblem& problem,
    const TransientHeat3DOptions& options
) {
    detail::validate_time_grid(options);
    const auto prepared = detail::prepare_transient_heat_3d(problem, options.start_time);
    detail::TransientHeat3DStepCache cache;
    return detail::solve_prepared_transient_heat_3d(
        prepared, options, prepared.initial_state, cache);
}

} // namespace nova::fem
