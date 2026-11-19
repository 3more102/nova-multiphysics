#include "nova/fem/transient_heat_3d.hpp"

#include "nova/fem/simplex_quadrature.hpp"
#include "nova/linalg/conjugate_gradient.hpp"
#include "p1_integration_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace nova::fem {
namespace {

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

} // namespace

TransientHeat3DResult solve_transient_heat_3d(
    const TransientHeat3DProblem& problem, const TransientHeat3DOptions& options
) {
    if (!std::isfinite(options.time_step) || !(options.time_step > 0.0) ||
        !std::isfinite(options.start_time) || options.steps == 0U) {
        throw std::invalid_argument("Transient 3D time step, start and step count are invalid.");
    }
    const double final_time = options.start_time +
        options.time_step * static_cast<double>(options.steps);
    if (!std::isfinite(final_time) || !(final_time > options.start_time)) {
        throw std::invalid_argument("Transient 3D final time must advance and be finite.");
    }
    if (options.scheme != TransientHeat3DScheme::BackwardEuler &&
        options.scheme != TransientHeat3DScheme::CrankNicolson) {
        throw std::invalid_argument("Unknown transient 3D time scheme.");
    }
    const double theta =
        options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 0.5;
    const double dt = options.time_step;

    const auto& diffusion = problem.diffusion;
    const auto& mesh = diffusion.mesh;
    // Validates geometry, connectivity, and exterior boundary assignments.
    detail::ExteriorFaceAssignments3D exterior(mesh);
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
    for (const double c : problem.cell_volumetric_capacity) {
        (void)positive(c, "Volumetric heat capacity");
    }
    for (const double q : diffusion.cell_source) {
        if (!std::isfinite(q)) throw std::invalid_argument("Transient 3D cell source must be finite.");
    }
    if (problem.initial_temperature.size() != n) {
        throw std::invalid_argument("Transient 3D initial temperature length mismatch.");
    }

    const detail::DiffusionMaterials3D materials(diffusion);
    std::vector<std::optional<double>> prescribed(n);
    std::set<std::size_t> explicit_nodes;
    for (const auto& bc : diffusion.dirichlet_constraints) {
        if (bc.node >= n) throw std::out_of_range("Transient 3D Dirichlet node is missing.");
        if (!explicit_nodes.insert(bc.node).second) {
            throw std::invalid_argument("Duplicate transient 3D Dirichlet vertex constraint.");
        }
        impose(prescribed, bc.node, bc.value);
    }
    for (const auto& bc : diffusion.dirichlet_faces) {
        exterior.claim(bc.face);
        if (!bc.value) throw std::invalid_argument("Transient 3D Dirichlet face needs a value.");
        for (const auto node : bc.face) {
            const auto& p = mesh.nodes[node];
            impose(prescribed, node, bc.value(p.x, p.y, p.z));
        }
    }
    for (const auto& bc : diffusion.neumann_faces) {
        exterior.claim(bc.face);
        if (!bc.outward_flux) {
            throw std::invalid_argument("Transient 3D Neumann face needs a flux.");
        }
    }
    for (const auto& bc : diffusion.robin_faces) {
        exterior.claim(bc.face);
        if (!bc.transfer_coefficient || !bc.ambient_value) {
            throw std::invalid_argument("Transient 3D Robin face needs transfer and ambient fields.");
        }
    }

    std::vector<double> state = problem.initial_temperature;
    for (std::size_t node = 0; node < n; ++node) {
        if (!std::isfinite(state[node])) {
            throw std::invalid_argument("Transient 3D initial temperatures must be finite.");
        }
        if (prescribed[node]) state[node] = *prescribed[node];
    }

    std::vector<Triplet> mass_entries, stiffness_entries, robin_entries;
    if (mesh.elements.size() > mass_entries.max_size() / 16U ||
        mesh.elements.size() > stiffness_entries.max_size() / 16U ||
        diffusion.robin_faces.size() > robin_entries.max_size() / 27U) {
        throw std::length_error("Transient 3D matrices exceed storage limits.");
    }
    mass_entries.reserve(mesh.elements.size() * 16U);
    stiffness_entries.reserve(mesh.elements.size() * 16U + diffusion.robin_faces.size() * 27U);
    robin_entries.reserve(diffusion.robin_faces.size() * 27U);
    std::vector<double> volume_load(n, 0.0), neumann_load(n, 0.0), ambient_load(n, 0.0);
    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        const auto& cell = mesh.elements[e];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, e);
        std::array<std::array<double, 4>, 4> local_mass{};
        std::array<double, 4> local_source{};
        const double extra_source = diffusion.cell_source.empty() ? 0.0 : diffusion.cell_source[e];
        const auto local_k = detail::integrate_p1_volume_3d(
            mesh, cell, basis, materials.on_cell(e),
            [&](const Point3D& p, const auto& shape, double weight) {
                const double capacity = cell_capacity
                    ? problem.cell_volumetric_capacity[e]
                    : positive(problem.volumetric_capacity(p.x, p.y, p.z),
                               "Volumetric heat capacity");
                const double spatial = diffusion.source
                    ? diffusion.source(p.x, p.y, p.z) : 0.0;
                if (!std::isfinite(spatial)) {
                    throw std::invalid_argument("Transient 3D source must be finite.");
                }
                const double source = spatial + extra_source;
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
            volume_load[cell[i]] += local_source[i];
            if (!std::isfinite(volume_load[cell[i]])) {
                throw std::overflow_error("Transient 3D source load overflowed.");
            }
            for (std::size_t j = 0; j < 4; ++j) {
                if (!std::isfinite(local_mass[i][j]) || !std::isfinite(local_k[i][j])) {
                    throw std::overflow_error("Transient 3D element matrix overflowed.");
                }
                mass_entries.push_back({cell[i], cell[j], local_mass[i][j]});
                stiffness_entries.push_back({cell[i], cell[j], local_k[i][j]});
            }
        }
    }

    for (const auto& bc : diffusion.neumann_faces) {
        const double area_scale = detail::face_jacobian(mesh, bc.face);
        for (const auto& q : triangle_face_degree_two_quadrature) {
            const auto p = detail::simplex_point(mesh, bc.face, q.shape);
            const double flux = bc.outward_flux(p.x, p.y, p.z);
            if (!std::isfinite(flux)) {
                throw std::invalid_argument("Transient 3D Neumann flux must be finite.");
            }
            for (std::size_t i = 0; i < 3; ++i) {
                neumann_load[bc.face[i]] -= flux * q.shape[i] * area_scale * q.weight;
                if (!std::isfinite(neumann_load[bc.face[i]])) {
                    throw std::overflow_error("Transient 3D Neumann load overflowed.");
                }
            }
        }
    }
    for (const auto& bc : diffusion.robin_faces) {
        const double area_scale = detail::face_jacobian(mesh, bc.face);
        for (const auto& q : triangle_face_degree_two_quadrature) {
            const auto p = detail::simplex_point(mesh, bc.face, q.shape);
            const double h = bc.transfer_coefficient(p.x, p.y, p.z);
            const double ambient = bc.ambient_value(p.x, p.y, p.z);
            if (!std::isfinite(h) || h < 0.0 || !std::isfinite(ambient)) {
                throw std::invalid_argument("Transient 3D Robin data must be finite with h >= 0.");
            }
            const double weight = area_scale * q.weight;
            for (std::size_t i = 0; i < 3; ++i) {
                ambient_load[bc.face[i]] += h * ambient * q.shape[i] * weight;
                if (!std::isfinite(ambient_load[bc.face[i]])) {
                    throw std::overflow_error("Transient 3D Robin load overflowed.");
                }
                for (std::size_t j = 0; j < 3; ++j) {
                    const double value = h * q.shape[i] * q.shape[j] * weight;
                    if (!std::isfinite(value)) {
                        throw std::overflow_error("Transient 3D Robin matrix overflowed.");
                    }
                    robin_entries.push_back({bc.face[i], bc.face[j], value});
                    stiffness_entries.push_back({bc.face[i], bc.face[j], value});
                }
            }
        }
    }

    std::vector<double> total_load(n);
    for (std::size_t i = 0; i < n; ++i) {
        total_load[i] = volume_load[i] + neumann_load[i] + ambient_load[i];
        if (!std::isfinite(total_load[i])) {
            throw std::overflow_error("Transient 3D total load overflowed.");
        }
    }
    const auto mass = CsrMatrix::from_triplets(n, n, mass_entries);
    const auto stiffness = CsrMatrix::from_triplets(n, n, stiffness_entries);
    const auto robin = CsrMatrix::from_triplets(n, n, robin_entries);
    const auto lhs = CsrMatrix::linear_combination(mass, 1.0, stiffness, theta * dt);
    const auto rhs_operator =
        CsrMatrix::linear_combination(mass, 1.0, stiffness, -(1.0 - theta) * dt);

    // Eliminate the constant Dirichlet block ONCE. The free/free CSR matrix
    // is reused for every step; the previous temperature warm-starts CG.
    std::vector<Triplet> free_entries;
    std::vector<Triplet> boundary_couplings;
    if (lhs.non_zeros() > free_entries.max_size() - n) {
        throw std::length_error("Transient 3D constrained matrix is too large.");
    }
    free_entries.reserve(lhs.non_zeros() + n);
    lhs.for_each_nonzero([&](std::size_t row, std::size_t col, double value) {
        if (!prescribed[row] && !prescribed[col]) {
            free_entries.push_back({row, col, value});
        } else if (!prescribed[row] && prescribed[col]) {
            boundary_couplings.push_back({row, col, value});
        }
    });
    for (std::size_t i = 0; i < n; ++i) {
        if (prescribed[i]) free_entries.push_back({i, i, 1.0});
    }
    const auto system = CsrMatrix::from_triplets(n, n, free_entries);
    std::vector<double> boundary_shift(n, 0.0);
    for (const auto& edge : boundary_couplings) {
        boundary_shift[edge.row] -= edge.value * *prescribed[edge.column];
    }
    const auto capacity_weights = mass.multiply(std::vector<double>(n, 1.0));
    TransientHeat3DResult result;
    result.final_time = final_time;
    result.linear_solvers.reserve(options.steps);
    result.diagnostics.reserve(options.steps);

    for (std::size_t step = 1; step <= options.steps; ++step) {
        std::vector<double> rhs = rhs_operator.multiply(state);
        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] += dt * total_load[i];
            if (!std::isfinite(rhs[i])) throw std::overflow_error("Transient 3D RHS overflowed.");
        }
        // Preserve the assembled, unconstrained rhs for reaction diagnostics.
        const auto physical_rhs = rhs;
        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] += boundary_shift[i];
            if (prescribed[i]) rhs[i] = *prescribed[i];
            if (!std::isfinite(rhs[i])) throw std::overflow_error("Transient 3D boundary RHS overflowed.");
        }
        // Solve for the small increment rather than the large absolute
        // temperature. A warm-started CG on absolute T defines its stopping
        // tolerance using a tiny initial residual, but true-residual
        // recomputation subtracts O(T) products and can hit the roundoff
        // floor. Defect correction keeps the same warm start while avoiding
        // that loss of relative precision inside CG.
        const auto current_product = system.multiply(state);
        std::vector<double> increment_rhs(n);
        for (std::size_t i = 0; i < n; ++i) {
            increment_rhs[i] = rhs[i] - current_product[i];
            if (!std::isfinite(increment_rhs[i])) {
                throw std::overflow_error("Transient 3D correction load overflowed.");
            }
        }
        auto solved = nova::linalg::solve_conjugate_gradient(
            system, increment_rhs, options.solver_options
        );
        if (!solved.converged) {
            throw std::runtime_error("Transient 3D step " + std::to_string(step) +
                                     " failed: " + solved.message);
        }
        std::vector<double> next(n);
        for (std::size_t i = 0; i < n; ++i) {
            next[i] = state[i] + solved.solution[i];
            if (prescribed[i]) next[i] = *prescribed[i];
            if (!std::isfinite(next[i])) {
                throw std::overflow_error("Transient 3D corrected temperature overflowed.");
            }
        }
        // Publish the physical state rather than the internal correction.
        solved.solution = next;
        const auto unconstrained_lhs = lhs.multiply(next);
        const auto previous_robin = robin.multiply(state);
        const auto next_robin = robin.multiply(next);

        TransientHeatStepDiagnostics3D diag;
        diag.time = options.start_time + dt * static_cast<double>(step);
        if (!std::isfinite(diag.time)) throw std::overflow_error("Transient 3D time overflowed.");
        diag.heat_content = dot(capacity_weights, next);
        diag.source_energy = dt * sum(volume_load);
        diag.neumann_energy = dt * sum(neumann_load);
        diag.robin_energy = dt * (
            sum(ambient_load) - theta * sum(next_robin) -
            (1.0 - theta) * sum(previous_robin)
        );
        for (std::size_t i = 0; i < n; ++i) {
            if (prescribed[i]) {
                diag.dirichlet_energy += unconstrained_lhs[i] - physical_rhs[i];
            }
        }
        diag.balance_residual = diag.heat_content - dot(capacity_weights, state) -
            diag.source_energy - diag.neumann_energy - diag.robin_energy -
            diag.dirichlet_energy;
        if (!std::isfinite(diag.dirichlet_energy) || !std::isfinite(diag.balance_residual)) {
            throw std::overflow_error("Transient 3D energy diagnostic overflowed.");
        }

        state = next;
        result.diagnostics.push_back(diag);
        result.linear_solvers.push_back(std::move(solved));
    }
    result.temperature = std::move(state);
    return result;
}

} // namespace nova::fem
