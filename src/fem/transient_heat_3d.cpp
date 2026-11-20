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
    const bool time_dependent_capacity =
        static_cast<bool>(problem.time_dependent_volumetric_capacity);
    if (time_dependent_capacity) {
        if (spatial_capacity || cell_capacity) {
            throw std::invalid_argument(
                "Time-dependent capacity cannot be combined with static capacity.");
        }
    } else if (spatial_capacity == cell_capacity) {
        throw std::invalid_argument("Choose exactly one transient 3D capacity representation.");
    }
    if (cell_capacity && problem.cell_volumetric_capacity.size() != mesh.elements.size()) {
        throw std::invalid_argument("Transient 3D cell capacity length mismatch.");
    }
    const bool time_dependent_scalar_conductivity =
        static_cast<bool>(problem.time_dependent_diffusivity);
    const bool time_dependent_tensor_conductivity =
        static_cast<bool>(problem.time_dependent_diffusivity_tensor);
    if (time_dependent_scalar_conductivity && time_dependent_tensor_conductivity) {
        throw std::invalid_argument(
            "Select time-dependent scalar or tensor conductivity, not both.");
    }
    const bool time_dependent_conductivity =
        time_dependent_scalar_conductivity || time_dependent_tensor_conductivity;
    const bool static_conductivity =
        static_cast<bool>(diffusion.diffusivity) ||
        static_cast<bool>(diffusion.diffusivity_tensor) ||
        !diffusion.cell_diffusivity.empty() ||
        !diffusion.cell_diffusivity_tensors.empty();
    if (time_dependent_conductivity && static_conductivity) {
        throw std::invalid_argument(
            "Time-dependent conductivity cannot be combined with static conductivity.");
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

    std::optional<DiffusionMaterials3D> static_materials;
    if (!time_dependent_conductivity) {
        static_materials.emplace(diffusion);
    }
    PreparedTransientHeat3D prepared;
    prepared.mesh = &mesh;
    prepared.time_dependent_source = problem.time_dependent_source;
    prepared.time_dependent_diffusivity = problem.time_dependent_diffusivity;
    prepared.time_dependent_diffusivity_tensor =
        problem.time_dependent_diffusivity_tensor;
    prepared.time_dependent_volumetric_capacity =
        problem.time_dependent_volumetric_capacity;
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
            const double jacobian = std::abs(basis.signed_jacobian_determinant);
            if (!std::isfinite(jacobian) || !(jacobian > 0.0)) {
                throw std::overflow_error("Transient 3D element volume is not representable.");
            }
            prepared.element_jacobians.push_back(jacobian);
        }
        std::array<std::array<double, 4>, 4> local_mass{};
        std::array<double, 4> local_source{};
        const double extra_source =
            diffusion.cell_source.empty() ? 0.0 : diffusion.cell_source[element];
        if (!time_dependent_conductivity) {
            // Keep the established static-material integration path unchanged.
            const auto local_stiffness = integrate_p1_volume_3d(
                mesh, cell, basis, static_materials->on_cell(element),
                [&](const Point3D& point, const auto& shape, double weight) {
                    if (!time_dependent_capacity) {
                        const double capacity = cell_capacity
                            ? problem.cell_volumetric_capacity[element]
                            : positive(problem.volumetric_capacity(
                                           point.x, point.y, point.z),
                                       "Volumetric heat capacity");
                        for (std::size_t i = 0; i < 4; ++i) {
                            for (std::size_t j = 0; j < 4; ++j) {
                                local_mass[i][j] += capacity * shape[i] * shape[j] * weight;
                            }
                        }
                    }
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
                    }
                });
            for (std::size_t i = 0; i < 4; ++i) {
                prepared.volume_load[cell[i]] += local_source[i];
                if (!std::isfinite(prepared.volume_load[cell[i]])) {
                    throw std::overflow_error("Transient 3D source load overflowed.");
                }
                for (std::size_t j = 0; j < 4; ++j) {
                    if (!std::isfinite(local_stiffness[i][j])) {
                        throw std::overflow_error("Transient 3D stiffness matrix overflowed.");
                    }
                    stiffness_entries.push_back(
                        {cell[i], cell[j], local_stiffness[i][j]});
                    if (!time_dependent_capacity) {
                        if (!std::isfinite(local_mass[i][j])) {
                            throw std::overflow_error("Transient 3D mass matrix overflowed.");
                        }
                        mass_entries.push_back({cell[i], cell[j], local_mass[i][j]});
                    }
                }
            }
        } else {
            // With dynamic conductivity only the static mass/source terms are
            // prepared here; K(x,y,z,t) is assembled at each time level below.
            const double jacobian = std::abs(basis.signed_jacobian_determinant);
            if (!std::isfinite(jacobian) || !(jacobian > 0.0)) {
                throw std::overflow_error("Transient 3D element volume is not representable.");
            }
            for (const auto& quadrature : tetrahedron_degree_two_quadrature) {
                const auto point = simplex_point(mesh, cell, quadrature.shape);
                const double weight = jacobian * quadrature.weight;
                if (!time_dependent_capacity) {
                    const double capacity = cell_capacity
                        ? problem.cell_volumetric_capacity[element]
                        : positive(problem.volumetric_capacity(
                                       point.x, point.y, point.z),
                                   "Volumetric heat capacity");
                    for (std::size_t i = 0; i < 4; ++i) {
                        for (std::size_t j = 0; j < 4; ++j) {
                            local_mass[i][j] +=
                                capacity * quadrature.shape[i] * quadrature.shape[j] * weight;
                        }
                    }
                }
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
                    local_source[i] += source * quadrature.shape[i] * weight;
                }
            }
            for (std::size_t i = 0; i < 4; ++i) {
                prepared.volume_load[cell[i]] += local_source[i];
                if (!std::isfinite(prepared.volume_load[cell[i]])) {
                    throw std::overflow_error("Transient 3D source load overflowed.");
                }
                if (!time_dependent_capacity) {
                    for (std::size_t j = 0; j < 4; ++j) {
                        if (!std::isfinite(local_mass[i][j])) {
                            throw std::overflow_error("Transient 3D mass matrix overflowed.");
                        }
                        mass_entries.push_back({cell[i], cell[j], local_mass[i][j]});
                    }
                }
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

void prepare_step_operators(
    const PreparedTransientHeat3D& prepared,
    double time_step,
    TransientHeat3DScheme scheme,
    const TransientHeat3DMaterialState& previous_material,
    const TransientHeat3DMaterialState& next_material,
    TransientHeat3DStepCache& cache,
    const CsrMatrix* previous_time_dependent_robin,
    const CsrMatrix* next_time_dependent_robin
) {
    const auto matches_cached_matrix = [](
        const std::optional<CsrMatrix>& cached,
        const CsrMatrix& candidate
    ) {
        return cached.has_value() && cached->exactly_equal(candidate);
    };
    if (cache.operators && cache.time_step == time_step &&
        cache.scheme == scheme &&
        matches_cached_matrix(cache.previous_material_mass, previous_material.mass) &&
        matches_cached_matrix(cache.next_material_mass, next_material.mass) &&
        matches_cached_matrix(cache.previous_material_stiffness, previous_material.stiffness) &&
        matches_cached_matrix(cache.next_material_stiffness, next_material.stiffness) &&
        (previous_time_dependent_robin == nullptr
            ? !cache.previous_time_dependent_robin.has_value()
            : matches_cached_matrix(cache.previous_time_dependent_robin,
                                    *previous_time_dependent_robin)) &&
        (next_time_dependent_robin == nullptr
            ? !cache.next_time_dependent_robin.has_value()
            : matches_cached_matrix(cache.next_time_dependent_robin,
                                    *next_time_dependent_robin))) {
        return;
    }

    const double theta =
        scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 0.5;
    const bool dynamic_conductivity =
        static_cast<bool>(prepared.time_dependent_diffusivity) ||
        static_cast<bool>(prepared.time_dependent_diffusivity_tensor);
    CsrMatrix lhs_stiffness = next_material.stiffness;
    CsrMatrix rhs_stiffness = previous_material.stiffness;
    if (dynamic_conductivity) {
        // Dynamic material states contain K(t) only; add static Robin terms
        // separately. The legacy static path already includes static Robin.
        lhs_stiffness = CsrMatrix::linear_combination(
            lhs_stiffness, 1.0, prepared.robin, 1.0);
        rhs_stiffness = CsrMatrix::linear_combination(
            rhs_stiffness, 1.0, prepared.robin, 1.0);
    }
    if (next_time_dependent_robin != nullptr) {
        lhs_stiffness = CsrMatrix::linear_combination(
            lhs_stiffness, 1.0, *next_time_dependent_robin, 1.0);
    }
    if (previous_time_dependent_robin != nullptr) {
        rhs_stiffness = CsrMatrix::linear_combination(
            rhs_stiffness, 1.0, *previous_time_dependent_robin, 1.0);
    }

    TransientHeat3DStepOperators next;
    next.lhs = CsrMatrix::linear_combination(
        next_material.mass, 1.0, lhs_stiffness, theta * time_step);
    next.rhs_operator = CsrMatrix::linear_combination(
        previous_material.mass, 1.0, rhs_stiffness, -(1.0 - theta) * time_step);

    const std::size_t n = prepared.initial_state.size();
    std::vector<Triplet> free_entries;
    std::vector<Triplet> boundary_couplings;
    if (n > free_entries.max_size() ||
        next.lhs.non_zeros() > free_entries.max_size() - n) {
        throw std::length_error("Transient 3D constrained matrix is too large.");
    }
    free_entries.reserve(next.lhs.non_zeros() + n);
    next.lhs.for_each_nonzero(
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
    next.system = CsrMatrix::from_triplets(n, n, free_entries);
    next.boundary_couplings = std::move(boundary_couplings);

    if (cache.operator_builds == std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error("Transient 3D operator-build counter overflowed.");
    }
    cache.operators = std::move(next);
    cache.time_step = time_step;
    cache.scheme = scheme;
    cache.previous_material_mass = previous_material.mass;
    cache.next_material_mass = next_material.mass;
    cache.previous_material_stiffness = previous_material.stiffness;
    cache.next_material_stiffness = next_material.stiffness;
    if (previous_time_dependent_robin != nullptr) {
        cache.previous_time_dependent_robin = *previous_time_dependent_robin;
    } else {
        cache.previous_time_dependent_robin.reset();
    }
    if (next_time_dependent_robin != nullptr) {
        cache.next_time_dependent_robin = *next_time_dependent_robin;
    } else {
        cache.next_time_dependent_robin.reset();
    }
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

struct TimeDependentRobinEvaluation3D {
    CsrMatrix matrix;
    std::vector<double> ambient_load;
};

TimeDependentRobinEvaluation3D evaluate_time_dependent_robin(
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
    std::vector<Triplet> entries;
    if (prepared.time_dependent_robin_faces.size() >
        entries.max_size() / 27U) {
        throw std::length_error(
            "Transient 3D time-dependent Robin matrix exceeds storage limits.");
    }
    entries.reserve(prepared.time_dependent_robin_faces.size() * 27U);
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
            const double weight = boundary.area_scale * quadrature.weight;
            const double ambient_transfer = transfer * ambient;
            if (!std::isfinite(weight) || !std::isfinite(ambient_transfer)) {
                throw std::overflow_error(
                    "Transient 3D time-dependent Robin quadrature overflowed.");
            }
            for (std::size_t i = 0; i < 3; ++i) {
                const auto row = boundary.face[i];
                const double load_contribution =
                    ambient_transfer * quadrature.shape[i] * weight;
                if (!std::isfinite(load_contribution)) {
                    throw std::overflow_error(
                        "Transient 3D time-dependent Robin load overflowed.");
                }
                ambient_load[row] += load_contribution;
                if (!std::isfinite(ambient_load[row])) {
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
                    entries.push_back({row, boundary.face[j], value});
                }
            }
        }
    }

    return {
        CsrMatrix::from_triplets(n, n, entries),
        std::move(ambient_load)
    };
}

TransientHeat3DMaterialState evaluate_time_dependent_materials(
    const PreparedTransientHeat3D& prepared,
    double time
) {
    if (!std::isfinite(time) || prepared.mesh == nullptr) {
        throw std::logic_error(
            "Transient 3D time-dependent material preparation is inconsistent.");
    }

    TransientHeat3DMaterialState state{
        prepared.mass, prepared.stiffness, prepared.capacity_weights
    };
    const bool dynamic_capacity =
        static_cast<bool>(prepared.time_dependent_volumetric_capacity);
    const bool dynamic_conductivity =
        static_cast<bool>(prepared.time_dependent_diffusivity) ||
        static_cast<bool>(prepared.time_dependent_diffusivity_tensor);
    if (!dynamic_capacity && !dynamic_conductivity) {
        return state;
    }

    SpatialScalarField3D scalar_field;
    SymmetricTensorField3D tensor_field;
    if (prepared.time_dependent_diffusivity) {
        const auto field = prepared.time_dependent_diffusivity;
        scalar_field = [field, time](double x, double y, double z) {
            return field(x, y, z, time);
        };
    }
    if (prepared.time_dependent_diffusivity_tensor) {
        const auto field = prepared.time_dependent_diffusivity_tensor;
        tensor_field = [field, time](double x, double y, double z) {
            return field(x, y, z, time);
        };
    }
    std::optional<ConductivityField3D> conductivity;
    if (dynamic_conductivity) {
        conductivity.emplace(ConductivityField3D{scalar_field, tensor_field});
        conductivity->validate();
    }

    const auto& mesh = *prepared.mesh;
    const std::size_t n = prepared.initial_state.size();
    std::vector<Triplet> mass_entries;
    std::vector<Triplet> stiffness_entries;
    if (dynamic_capacity && mesh.elements.size() > mass_entries.max_size() / 16U) {
        throw std::length_error(
            "Transient 3D time-dependent mass matrix exceeds storage limits.");
    }
    if (dynamic_conductivity &&
        mesh.elements.size() > stiffness_entries.max_size() / 16U) {
        throw std::length_error(
            "Transient 3D time-dependent stiffness matrix exceeds storage limits.");
    }
    if (dynamic_capacity) mass_entries.reserve(mesh.elements.size() * 16U);
    if (dynamic_conductivity) stiffness_entries.reserve(mesh.elements.size() * 16U);

    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const auto& cell = mesh.elements[element];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, element);
        std::array<std::array<double, 4>, 4> local_mass{};
        P1Matrix3D local_stiffness{};
        if (conductivity) {
            local_stiffness = integrate_p1_volume_3d(
                mesh, cell, basis, *conductivity,
                [&](const Point3D& point, const auto& shape, double weight) {
                    if (dynamic_capacity) {
                        const double capacity = positive(
                            prepared.time_dependent_volumetric_capacity(
                                point.x, point.y, point.z, time),
                            "Time-dependent volumetric heat capacity");
                        for (std::size_t i = 0; i < 4; ++i) {
                            for (std::size_t j = 0; j < 4; ++j) {
                                local_mass[i][j] += capacity * shape[i] * shape[j] * weight;
                            }
                        }
                    }
                });
        } else if (dynamic_capacity) {
            const double jacobian = std::abs(basis.signed_jacobian_determinant);
            if (!std::isfinite(jacobian) || !(jacobian > 0.0)) {
                throw std::overflow_error(
                    "Transient 3D material element volume is not representable.");
            }
            for (const auto& quadrature : tetrahedron_degree_two_quadrature) {
                const auto point = simplex_point(mesh, cell, quadrature.shape);
                const double capacity = positive(
                    prepared.time_dependent_volumetric_capacity(
                        point.x, point.y, point.z, time),
                    "Time-dependent volumetric heat capacity");
                const double weight = jacobian * quadrature.weight;
                for (std::size_t i = 0; i < 4; ++i) {
                    for (std::size_t j = 0; j < 4; ++j) {
                        local_mass[i][j] +=
                            capacity * quadrature.shape[i] * quadrature.shape[j] * weight;
                    }
                }
            }
        }
        if (dynamic_capacity) {
            for (std::size_t i = 0; i < 4; ++i) {
                for (std::size_t j = 0; j < i; ++j) {
                    local_mass[j][i] = local_mass[i][j];
                }
            }
        }
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = 0; j < 4; ++j) {
                if (dynamic_capacity) {
                    if (!std::isfinite(local_mass[i][j])) {
                        throw std::overflow_error(
                            "Transient 3D time-dependent mass matrix overflowed.");
                    }
                    mass_entries.push_back({cell[i], cell[j], local_mass[i][j]});
                }
                if (dynamic_conductivity) {
                    if (!std::isfinite(local_stiffness[i][j])) {
                        throw std::overflow_error(
                            "Transient 3D time-dependent stiffness matrix overflowed.");
                    }
                    stiffness_entries.push_back(
                        {cell[i], cell[j], local_stiffness[i][j]});
                }
            }
        }
    }

    if (dynamic_capacity) {
        state.mass = CsrMatrix::from_triplets(n, n, mass_entries);
        state.capacity_weights = state.mass.multiply(
            std::vector<double>(n, 1.0));
    }
    if (dynamic_conductivity) {
        state.stiffness = CsrMatrix::from_triplets(n, n, stiffness_entries);
    }
    return state;
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

    const double dt = options.time_step;
    const double theta =
        options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 0.5;
    auto previous_material =
        evaluate_time_dependent_materials(prepared, options.start_time);
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
    std::optional<TimeDependentRobinEvaluation3D> previous_time_robin;
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
        std::optional<TimeDependentRobinEvaluation3D> next_time_robin;
        if (!prepared.time_dependent_robin_faces.empty()) {
            next_time_robin = evaluate_time_dependent_robin(prepared, step_time);
        }
        auto next_material =
            evaluate_time_dependent_materials(prepared, step_time);
        prepare_step_operators(
            prepared, dt, options.scheme, previous_material, next_material, cache,
            previous_time_robin ? &previous_time_robin->matrix : nullptr,
            next_time_robin ? &next_time_robin->matrix : nullptr);
        const auto& operators = *cache.operators;

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
        std::vector<double> weighted_time_robin_load;
        if (!prepared.time_dependent_robin_faces.empty()) {
            if (!next_time_robin || (theta < 1.0 && !previous_time_robin)) {
                throw std::logic_error(
                    "Transient 3D time-dependent Robin endpoint data are incomplete.");
            }
            weighted_time_robin_load.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                weighted_time_robin_load[i] =
                    theta * next_time_robin->ambient_load[i];
                if (theta < 1.0) {
                    weighted_time_robin_load[i] +=
                        (1.0 - theta) * previous_time_robin->ambient_load[i];
                }
                if (!std::isfinite(weighted_time_robin_load[i])) {
                    throw std::overflow_error(
                        "Transient 3D weighted Robin load overflowed.");
                }
            }
        }

        const auto current_prescribed =
            evaluate_prescribed_values(prepared, step_time);

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
        std::vector<double> rhs = operators.rhs_operator.multiply(state);
        for (std::size_t i = 0; i < n; ++i) {
            const double load = prepared.total_load[i] +
                (prepared.time_dependent_source ? weighted_time_source[i] : 0.0) +
                (!prepared.time_dependent_neumann_faces.empty()
                    ? weighted_time_neumann[i] : 0.0) +
                (!prepared.time_dependent_robin_faces.empty()
                    ? weighted_time_robin_load[i] : 0.0);
            rhs[i] += dt * load;
            if (!std::isfinite(rhs[i])) {
                throw std::overflow_error("Transient 3D RHS overflowed.");
            }
        }
        // Preserve the assembled, unconstrained rhs for reaction diagnostics.
        const auto physical_rhs = rhs;
        std::vector<double> boundary_shift(n, 0.0);
        for (const auto& edge : operators.boundary_couplings) {
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
        const auto current_product = operators.system.multiply(state);
        std::vector<double> increment_rhs(n);
        for (std::size_t i = 0; i < n; ++i) {
            increment_rhs[i] = rhs[i] - current_product[i];
            if (!std::isfinite(increment_rhs[i])) {
                throw std::overflow_error("Transient 3D correction load overflowed.");
            }
        }
        auto solved = cache.cached_cg.solve(
            operators.system, increment_rhs, options.solver_options);
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

        const auto unconstrained_lhs = operators.lhs.multiply(next);
        const auto previous_robin = prepared.robin.multiply(state);
        const auto next_robin = prepared.robin.multiply(next);
        const auto previous_time_robin_product = previous_time_robin
            ? previous_time_robin->matrix.multiply(state)
            : std::vector<double>{};
        const auto next_time_robin_product = next_time_robin
            ? next_time_robin->matrix.multiply(next)
            : std::vector<double>{};

        TransientHeatStepDiagnostics3D diagnostic;
        diagnostic.time = step_time;
        if (!std::isfinite(diagnostic.time)) {
            throw std::overflow_error("Transient 3D time overflowed.");
        }
        diagnostic.heat_content = dot(next_material.capacity_weights, next);
        diagnostic.source_energy = dt * (
            sum(prepared.volume_load) + sum(weighted_time_source));
        if (!std::isfinite(diagnostic.source_energy)) {
            throw std::overflow_error("Transient 3D source energy overflowed.");
        }
        diagnostic.neumann_energy = dt * (
            sum(prepared.neumann_load) + sum(weighted_time_neumann));
        diagnostic.robin_energy = dt * (
            sum(prepared.ambient_load) + sum(weighted_time_robin_load) -
            theta * (sum(next_robin) + sum(next_time_robin_product)) -
            (1.0 - theta) *
                (sum(previous_robin) + sum(previous_time_robin_product)));
        for (std::size_t i = 0; i < n; ++i) {
            if (current_prescribed[i]) {
                diagnostic.dirichlet_energy +=
                    unconstrained_lhs[i] - physical_rhs[i];
            }
        }
        diagnostic.balance_residual =
            diagnostic.heat_content - dot(previous_material.capacity_weights, state) -
            diagnostic.source_energy - diagnostic.neumann_energy -
            diagnostic.robin_energy - diagnostic.dirichlet_energy;
        if (!std::isfinite(diagnostic.source_energy) ||
            !std::isfinite(diagnostic.neumann_energy) ||
            !std::isfinite(diagnostic.robin_energy) ||
            !std::isfinite(diagnostic.dirichlet_energy) ||
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
        previous_material = std::move(next_material);
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
