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

double positive(double value, const char* what) {
    if (!std::isfinite(value) || !(value > 0.0)) {
        throw std::invalid_argument(std::string(what) + " must be finite and positive.");
    }
    return value;
}

void impose(std::vector<std::optional<double>>& constraints,
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

double sum(const std::vector<double>& values) {
    long double total = 0.0L;
    for (double value : values) total += value;
    if (!std::isfinite(total) ||
        std::abs(total) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Transient 3D energy sum is outside the finite range.");
    }
    return static_cast<double>(total);
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
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

void validate_time_grid(const TransientHeat3DOptions& options) {
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
    const TransientHeat3DProblem& problem
) {
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
    prepared.prescribed.resize(n);
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

    prepared.initial_state = problem.initial_temperature;
    for (std::size_t node = 0; node < n; ++node) {
        if (!std::isfinite(prepared.initial_state[node])) {
            throw std::invalid_argument(
                "Transient 3D initial temperatures must be finite.");
        }
        if (prepared.prescribed[node]) {
            prepared.initial_state[node] = *prepared.prescribed[node];
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

    const std::size_t n = prepared.initial_state.size();
    std::vector<Triplet> free_entries;
    std::vector<Triplet> boundary_couplings;
    if (n > free_entries.max_size() ||
        next.lhs.non_zeros() > free_entries.max_size() - n ||
        next.lhs.non_zeros() > boundary_couplings.max_size()) {
        throw std::length_error("Transient 3D constrained matrix is too large.");
    }
    free_entries.reserve(next.lhs.non_zeros() + n);
    boundary_couplings.reserve(next.lhs.non_zeros());
    next.lhs.for_each_nonzero(
        [&](std::size_t row, std::size_t column, double value) {
            if (!prepared.prescribed[row] && !prepared.prescribed[column]) {
                free_entries.push_back({row, column, value});
            } else if (!prepared.prescribed[row] &&
                       prepared.prescribed[column]) {
                boundary_couplings.push_back({row, column, value});
            }
        });
    for (std::size_t i = 0; i < n; ++i) {
        if (prepared.prescribed[i]) {
            free_entries.push_back({i, i, 1.0});
        }
    }
    next.system = CsrMatrix::from_triplets(n, n, free_entries);
    next.boundary_shift.assign(n, 0.0);
    for (const auto& edge : boundary_couplings) {
        next.boundary_shift[edge.row] -=
            edge.value * *prepared.prescribed[edge.column];
        if (!std::isfinite(next.boundary_shift[edge.row])) {
            throw std::overflow_error("Transient 3D boundary shift overflowed.");
        }
    }

    if (cache.operator_builds == std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error("Transient 3D operator-build counter overflowed.");
    }
    cache.operators = std::move(next);
    cache.time_step = time_step;
    cache.scheme = scheme;
    ++cache.operator_builds;
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
    for (std::size_t node = 0; node < n; ++node) {
        if (!std::isfinite(state[node])) {
            throw std::invalid_argument("Transient 3D initial temperatures must be finite.");
        }
        if (prepared.prescribed[node]) {
            state[node] = *prepared.prescribed[node];
        }
    }

    prepare_step_operators(prepared, options.time_step, options.scheme, cache);
    const auto& operators = *cache.operators;
    const double dt = options.time_step;
    const double theta =
        options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 0.5;

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
        std::vector<double> rhs = operators.rhs_operator.multiply(state);
        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] += dt * prepared.total_load[i];
            if (!std::isfinite(rhs[i])) {
                throw std::overflow_error("Transient 3D RHS overflowed.");
            }
        }
        // Preserve the assembled, unconstrained rhs for reaction diagnostics.
        const auto physical_rhs = rhs;
        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] += operators.boundary_shift[i];
            if (prepared.prescribed[i]) {
                rhs[i] = *prepared.prescribed[i];
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
            if (prepared.prescribed[i]) {
                next[i] = *prepared.prescribed[i];
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

        TransientHeatStepDiagnostics3D diagnostic;
        diagnostic.time = options.start_time +
            dt * static_cast<double>(step);
        if (!std::isfinite(diagnostic.time)) {
            throw std::overflow_error("Transient 3D time overflowed.");
        }
        diagnostic.heat_content = dot(prepared.capacity_weights, next);
        diagnostic.source_energy = dt * sum(prepared.volume_load);
        diagnostic.neumann_energy = dt * sum(prepared.neumann_load);
        diagnostic.robin_energy = dt * (
            sum(prepared.ambient_load) - theta * sum(next_robin) -
            (1.0 - theta) * sum(previous_robin));
        for (std::size_t i = 0; i < n; ++i) {
            if (prepared.prescribed[i]) {
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
    const auto prepared = detail::prepare_transient_heat_3d(problem);
    detail::TransientHeat3DStepCache cache;
    return detail::solve_prepared_transient_heat_3d(
        prepared, options, prepared.initial_state, cache);
}

} // namespace nova::fem
