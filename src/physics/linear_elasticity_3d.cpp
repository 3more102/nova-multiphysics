#include "nova/physics/thermoelasticity_3d.hpp"

#include "nova/fem/tetrahedron_h1.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::physics {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed(Clock::time_point first, Clock::time_point last) {
    return std::chrono::duration<double>(last - first).count();
}

struct LameParameters {
    double lambda{};
    double mu{};
    double bulk{}; // K = E / (3(1-2 nu))
};

[[nodiscard]] LameParameters lame(const IsotropicElasticMaterial3D& material) {
    const double e = material.young_modulus;
    const double nu = material.poisson_ratio;
    if (!std::isfinite(e) || !(e > 0.0) ||
        !std::isfinite(nu) || !(nu > -1.0 && nu < 0.5)) {
        throw std::invalid_argument(
            "3D isotropic elasticity requires finite E > 0 and -1 < nu < 0.5.");
    }
    const double mu = e / (2.0 * (1.0 + nu));
    const double lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double bulk = e / (3.0 * (1.0 - 2.0 * nu));
    if (!std::isfinite(mu) || !(mu > 0.0) || !std::isfinite(lambda) ||
        !std::isfinite(bulk) || !(bulk > 0.0)) {
        throw std::overflow_error("3D elasticity Lame coefficients are not representable.");
    }
    return {lambda, mu, bulk};
}

[[nodiscard]] bool finite(const ElasticVector3D& x) {
    return std::all_of(x.begin(), x.end(),
                       [](double value) { return std::isfinite(value); });
}

[[nodiscard]] double norm_squared(const ElasticVector3D& a) {
    return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
}

[[nodiscard]] ElasticVector3D between(
    const fem::Point3D& a, const fem::Point3D& b
) {
    return {b.x - a.x, b.y - a.y, b.z - a.z};
}

[[nodiscard]] ElasticVector3D cross(
    const ElasticVector3D& a, const ElasticVector3D& b
) {
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

[[nodiscard]] fem::BoundaryFace3D sorted_face(fem::BoundaryFace3D face) {
    std::sort(face.begin(), face.end());
    return face;
}

using FixedComponents = std::array<std::optional<double>, 3>;
using RigidRow = std::array<double, 6>;

// For u=t+omega cross r, each prescribed scalar displacement contributes
// one row on the six rigid modes. Centering and scaling the anchor coordinates
// makes the rank test invariant to the model origin and length unit.
std::size_t rigid_constraint_rank(const fem::TetrahedralMesh3D& mesh,
                                  const std::vector<std::size_t>& nodes,
                                  std::span<const FixedComponents> fixed) {
    if (nodes.empty()) return 0;
    const auto& origin = mesh.nodes[nodes.front()];
    double scale = 0.0;
    for (const auto node : nodes) {
        const auto r = between(origin, mesh.nodes[node]);
        scale = std::max({scale, std::abs(r[0]), std::abs(r[1]), std::abs(r[2])});
    }
    if (scale == 0.0) scale = 1.0;
    std::array<RigidRow, 6> basis{};
    std::size_t rank = 0;
    for (const auto node : nodes) {
        auto r = between(origin, mesh.nodes[node]);
        for (double& v : r) v /= scale;
        const auto [x,y,z] = r;
        const std::array<RigidRow, 3> rows{{
            {1,0,0,0,z,-y}, {0,1,0,-z,0,x}, {0,0,1,y,-x,0}
        }};
        for (std::size_t i = 0; i < 3; ++i) {
            if (!fixed[node][i]) continue;
            auto row = rows[i];
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t k = 0; k < rank; ++k) {
                    const double dot = std::inner_product(row.begin(), row.end(), basis[k].begin(), 0.0);
                    for (std::size_t j = 0; j < 6; ++j) row[j] -= dot * basis[k][j];
                }
            }
            const double norm = std::sqrt(std::inner_product(row.begin(), row.end(), row.begin(), 0.0));
            if (norm <= 1e-12) continue;
            for (std::size_t j = 0; j < 6; ++j) basis[rank][j] = row[j]/norm;
            if (++rank == 6) return rank;
        }
    }
    return rank;
}

void require_rigid_mode_anchors(const fem::TetrahedralMesh3D& mesh,
                                std::span<const std::size_t> components,
                                std::span<const FixedComponents> fixed) {
    std::map<std::size_t, std::vector<std::size_t>> anchors;
    for (std::size_t e = 0; e < components.size(); ++e) {
        auto& nodes = anchors[components[e]];
        for (const auto node : mesh.elements[e]) {
            if (fixed[node][0] || fixed[node][1] || fixed[node][2]) nodes.push_back(node);
        }
    }
    for (auto& [body, nodes] : anchors) {
        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        const auto rank = rigid_constraint_rank(mesh, nodes, fixed);
        if (rank < 6) {
            throw std::invalid_argument("Elastic supports on the solid containing tetrahedron " +
                std::to_string(body) + " constrain " + std::to_string(rank) + " of 6 rigid modes.");
        }
    }
}

[[nodiscard]] ElasticTensor6 strain_from_gradient(
    const std::array<ElasticVector3D, 3>& grad
) {
    return {
        grad[0][0], grad[1][1], grad[2][2],
        0.5 * (grad[0][1] + grad[1][0]),
        0.5 * (grad[1][2] + grad[2][1]),
        0.5 * (grad[0][2] + grad[2][0])
    };
}

[[nodiscard]] ElasticTensor6 stress_from_strain(
    const ElasticTensor6& eps, const LameParameters& mat
) {
    const double mean = (eps[0] + eps[1] + eps[2]) / 3.0;
    const double volumetric = (mat.bulk * mean) * 3.0;
    return {
        volumetric + 2.0 * mat.mu * (eps[0] - mean),
        volumetric + 2.0 * mat.mu * (eps[1] - mean),
        volumetric + 2.0 * mat.mu * (eps[2] - mean),
        2.0 * mat.mu * eps[3],
        2.0 * mat.mu * eps[4],
        2.0 * mat.mu * eps[5]
    };
}

struct ElasticCell {
    fem::Tetrahedron3D nodes;
    std::array<ElasticVector3D, 4> gradients;
    double volume;
    LameParameters material;
};

struct CellMoment {
    double mean{};
    long double variance{};
};

using Block = std::array<ElasticVector3D, 3>;

void add_stiffness_block(Block& block, const ElasticCell& cell,
                         const ElasticVector3D& a, const ElasticVector3D& b) {
    const auto& m = cell.material;
    const double dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            block[i][j] += cell.volume * (m.lambda * a[i] * b[j] +
                m.mu * a[j] * b[i] + (i == j ? m.mu * dot : 0.0));
        }
    }
}

// Gather one block row at a time. Valence-sized scratch replaces the 144
// scalar triplets per tetrahedron and the second, map-based assembly pass.
linalg::CsrMatrix assemble_stiffness(
    std::span<const ElasticCell> cells,
    std::span<const FixedComponents> fixed,
    std::vector<double>& base_load
) {
    struct Incident { std::size_t element, local_node; };
    std::vector<std::vector<Incident>> incidents(fixed.size());
    for (std::size_t e = 0; e < cells.size(); ++e) {
        for (std::size_t a = 0; a < 4; ++a) {
            incidents[cells[e].nodes[a]].push_back({e, a});
        }
    }
    std::vector<std::size_t> offsets{0}, columns;
    std::vector<double> values;
    offsets.reserve(base_load.size() + 1U);
    std::vector<std::size_t> neighbors, slot(fixed.size());
    std::vector<Block> blocks;
    for (std::size_t node = 0; node < fixed.size(); ++node) {
        const bool all_fixed = fixed[node][0] && fixed[node][1] && fixed[node][2];
        neighbors.clear();
        if (!all_fixed) {
            for (const auto& inc : incidents[node]) {
                const auto& nodes = cells[inc.element].nodes;
                neighbors.insert(neighbors.end(), nodes.begin(), nodes.end());
            }
            std::sort(neighbors.begin(), neighbors.end());
            neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
            for (std::size_t k = 0; k < neighbors.size(); ++k) slot[neighbors[k]] = k;
            blocks.assign(neighbors.size(), {});
            for (const auto& inc : incidents[node]) {
                const auto& cell = cells[inc.element];
                const auto& a = cell.gradients[inc.local_node];
                for (std::size_t local = 0; local < 4; ++local) {
                    add_stiffness_block(blocks[slot[cell.nodes[local]]], cell, a, cell.gradients[local]);
                }
            }
        }
        for (std::size_t i = 0; i < 3; ++i) {
            const auto row = 3U * node + i;
            auto& load = base_load[row];
            if (fixed[node][i]) {
                columns.push_back(row);
                values.push_back(1.0);
                offsets.push_back(values.size());
                load = *fixed[node][i];
                continue;
            }
            for (std::size_t k = 0; k < neighbors.size(); ++k) {
                const auto other = neighbors[k];
                for (std::size_t j = 0; j < 3; ++j) {
                    const double value = blocks[k][i][j];
                    if (!std::isfinite(value)) {
                        throw std::overflow_error("Non-finite 3D elastic stiffness.");
                    }
                    if (fixed[other][j]) load -= value * *fixed[other][j];
                    else if (value != 0.0) {
                        columns.push_back(3U * other + j);
                        values.push_back(value);
                    }
                }
            }
            if (!std::isfinite(load)) throw std::overflow_error("Non-finite elastic load.");
            offsets.push_back(values.size());
        }
    }
    return linalg::CsrMatrix::from_sorted_rows(base_load.size(), base_load.size(),
        std::move(offsets), std::move(columns), std::move(values));
}

} // namespace

struct PreparedElasticity3D::Impl {
    std::vector<ElasticCell> cells;
    std::vector<FixedComponents> fixed;
    linalg::CsrMatrix stiffness;
    std::vector<double> base_load, guess;
    linalg::CachedConjugateGradient cg;
    bool has_converged_guess{};
    ElasticityTimings3D profile;

    LinearElasticity3DResult solve(std::span<const CellMoment> thermal,
                                   const linalg::SolverOptions& options, Clock::time_point load_start);
};

PreparedElasticity3D::PreparedElasticity3D(
    const fem::TetrahedralMesh3D& mesh, std::span<const IsotropicElasticMaterial3D> materials,
    std::span<const DisplacementConstraint3D> constraints,
    std::span<const ElasticVector3D> body, std::span<const ElasticFaceTraction3D> tractions,
    std::span<const DisplacementComponent3D> component_constraints
) : impl_(std::make_unique<Impl>()) {
    const auto start = Clock::now();
    const auto topology = fem::analyze_tetrahedral_topology_3d(mesh);
    const auto topology_done = Clock::now();
    if (mesh.nodes.size() > (std::numeric_limits<std::size_t>::max() - 1U) / 3U) {
        throw std::length_error("3D elasticity DOF count overflows.");
    }
    if (materials.size() != mesh.elements.size() ||
        (!body.empty() && body.size() != materials.size())) {
        throw std::invalid_argument("3D elasticity requires one material and optional body force per cell.");
    }
    auto& state = *impl_;
    state.fixed.resize(mesh.nodes.size());
    state.base_load.assign(3U * mesh.nodes.size(), 0.0);
    state.guess = state.base_load;
    const auto impose = [&](std::size_t node, std::size_t axis, double value) {
        if (node >= state.fixed.size() || axis >= 3) throw std::out_of_range("Invalid elastic constraint DOF.");
        if (state.fixed[node][axis] || !std::isfinite(value)) {
            throw std::invalid_argument("Duplicate prescribed elastic DOF or non-finite displacement.");
        }
        state.fixed[node][axis] = value;
        state.guess[3U * node + axis] = value;
    };
    for (const auto& bc : constraints) {
        for (std::size_t i = 0; i < 3; ++i) impose(bc.node, i, bc.value[i]);
    }
    for (const auto& bc : component_constraints) impose(bc.node, bc.component, bc.value);
    require_rigid_mode_anchors(mesh, topology.cell_component, state.fixed);
    state.cells.reserve(mesh.elements.size());
    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        const auto geometry = fem::evaluate_tetrahedron_p1_physical_gradients_3d(mesh, e);
        const auto volume = std::abs(geometry.signed_jacobian_determinant) / 6.0;
        state.cells.push_back({mesh.elements[e], geometry.gradients, volume, lame(materials[e])});
        if (!body.empty()) {
            if (!finite(body[e])) throw std::invalid_argument("Non-finite 3D body force.");
            for (const auto node : mesh.elements[e]) {
                for (std::size_t i = 0; i < 3; ++i) {
                    state.base_load[3U*node+i] += volume * body[e][i] / 4.0;
                }
            }
        }
    }
    if (!tractions.empty()) {
        std::set<fem::BoundaryFace3D> available;
        for (const auto& face : topology.exterior_faces) available.insert(sorted_face(face));
        for (const auto& bc : tractions) {
            if (!available.erase(sorted_face(bc.face)) || !finite(bc.traction)) {
                throw std::invalid_argument("Elastic traction requires a finite load and a distinct exterior face.");
            }
            const auto& p = mesh.nodes[bc.face[0]];
            const auto area = 0.5 * std::sqrt(norm_squared(cross(
                between(p, mesh.nodes[bc.face[1]]), between(p, mesh.nodes[bc.face[2]]))));
            if (!(area > 0.0) || !std::isfinite(area)) {
                throw std::overflow_error("Invalid 3D elastic traction face area.");
            }
            for (const auto node : bc.face) {
                for (std::size_t i = 0; i < 3; ++i) state.base_load[3U*node+i] += area * bc.traction[i] / 3.0;
            }
        }
    }
    for (const auto load : state.base_load) {
        if (!std::isfinite(load)) throw std::overflow_error("Non-finite 3D elastic applied load.");
    }
    const auto elements_done = Clock::now();
    state.stiffness = assemble_stiffness(state.cells, state.fixed, state.base_load);
    state.profile.topology = elapsed(start, topology_done);
    state.profile.element_setup = elapsed(topology_done, elements_done);
    state.profile.stiffness_assembly = elapsed(elements_done, Clock::now());
}

PreparedElasticity3D::~PreparedElasticity3D() = default;
PreparedElasticity3D::PreparedElasticity3D(PreparedElasticity3D&&) noexcept = default;
PreparedElasticity3D& PreparedElasticity3D::operator=(PreparedElasticity3D&&) noexcept = default;

ElasticityTimings3D PreparedElasticity3D::timings() const noexcept {
    return impl_->profile;
}

std::size_t PreparedElasticity3D::stiffness_nonzeros() const noexcept {
    return impl_->stiffness.non_zeros();
}

std::size_t PreparedElasticity3D::factor_builds() const noexcept {
    return impl_->cg.factor_builds();
}

LinearElasticity3DResult PreparedElasticity3D::Impl::solve(
    std::span<const CellMoment> thermal, const linalg::SolverOptions& options, Clock::time_point load_start
) {
    auto rhs = base_load;
    for (std::size_t e = 0; e < thermal.size(); ++e) {
        const auto& cell = cells[e];
        const double pressure = (cell.material.bulk * thermal[e].mean) * 3.0;
        if (!std::isfinite(pressure)) throw std::overflow_error("Non-finite thermal eigenstress.");
        for (std::size_t a = 0; a < 4; ++a) {
            const auto node = cell.nodes[a];
            for (std::size_t i = 0; i < 3; ++i) {
                if (fixed[node][i]) continue;
                double& load = rhs[3U * node + i];
                load += cell.volume * pressure * cell.gradients[a][i];
                if (!std::isfinite(load)) throw std::overflow_error("Non-finite thermal expansion load.");
            }
        }
    }
    LinearElasticity3DResult result;
    auto settings = options;
    std::span<const double> initial = guess;
    if (has_converged_guess && options.relative_tolerance > 0 &&
        std::isfinite(options.relative_tolerance) && options.absolute_tolerance >= 0 &&
        std::isfinite(options.absolute_tolerance)) {
        // Warm starts should meet the same free-force residual target as a
        // fresh solve. Scaling by the much smaller *previous* defect can
        // otherwise demand precision below roundoff at an unchanged load.
        double relative_target = 0;
        bool homogeneous = true;
        for (std::size_t node = 0; node < fixed.size(); ++node) {
            for (std::size_t i = 0; i < 3; ++i) {
                if (fixed[node][i]) continue;
                const double f = rhs[3U*node+i];
                homogeneous = homogeneous && (f == 0);
                relative_target = std::hypot(relative_target, options.relative_tolerance*f);
            }
        }
        if (homogeneous) initial = rhs; // Identity rows contain the prescribed displacements.
        settings.absolute_tolerance = std::max(options.absolute_tolerance, relative_target);
        if (settings.absolute_tolerance > 0) settings.relative_tolerance = 0;
    }
    const auto solve_start = Clock::now();
    result.linear_solver = cg.solve(stiffness, rhs, settings, initial);
    const auto solve_done = Clock::now();
    result.displacement.resize(fixed.size());
    for (std::size_t n = 0; n < fixed.size(); ++n) {
        for (std::size_t i = 0; i < 3; ++i) {
            result.displacement[n][i] = result.linear_solver.solution[3U*n+i];
        }
    }
    result.cell_strain.reserve(cells.size());
    result.cell_stress.reserve(cells.size());
    long double total_energy = 0;
    for (std::size_t e = 0; e < cells.size(); ++e) {
        const auto& cell = cells[e];
        std::array<ElasticVector3D, 3> gradient{};
        for (std::size_t a = 0; a < 4U; ++a) {
            for (std::size_t i = 0; i < 3U; ++i) {
                for (std::size_t j = 0; j < 3U; ++j) {
                    gradient[i][j] += result.displacement[cell.nodes[a]][i] * cell.gradients[a][j];
                }
            }
        }
        const auto strain = strain_from_gradient(gradient);
        const auto moment = thermal.empty() ? CellMoment{} : thermal[e];
        auto elastic = strain;
        for (std::size_t i = 0; i < 3; ++i) elastic[i] -= moment.mean;
        const auto stress = stress_from_strain(elastic, cell.material);
        // For affine thermal strain s(x), mean[(eps-sI):C:(eps-sI)]
        // equals the centroid work plus (I:C:I)*Var(s).
        long double work = 9.0L * moment.variance * cell.material.bulk;
        for (std::size_t i = 0; i < 6; ++i) {
            work += (i < 3 ? 1.0L : 2.0L) * elastic[i] * stress[i];
        }
        total_energy += 0.5L * cell.volume * work;
        result.cell_strain.push_back(strain);
        result.cell_stress.push_back(stress);
    }
    result.strain_energy_joules = static_cast<double>(total_energy);
    if (!std::isfinite(result.strain_energy_joules)) {
        throw std::overflow_error("Non-finite 3D elastic strain energy.");
    }
    if (result.linear_solver.converged) {
        guess = result.linear_solver.solution;
        has_converged_guess = true;
    }
    profile.thermal_load += elapsed(load_start, solve_start);
    profile.linear_solve += elapsed(solve_start, solve_done);
    profile.recovery += elapsed(solve_done, Clock::now());
    return result;
}

LinearElasticity3DResult PreparedElasticity3D::solve_temperature(
    std::span<const double> temperature, std::span<const ThermalExpansion3D> expansion,
    const linalg::SolverOptions& options
) {
    const auto load_start = Clock::now();
    if (temperature.size() != impl_->fixed.size() || expansion.size() != impl_->cells.size()) {
        throw std::invalid_argument("Thermoelasticity requires nodal temperatures and one expansion law per cell.");
    }
    std::vector<CellMoment> moments(expansion.size());
    for (std::size_t e = 0; e < expansion.size(); ++e) {
        const auto [alpha, reference] = expansion[e];
        if (!std::isfinite(alpha) || !std::isfinite(reference)) {
            throw std::invalid_argument("Thermal expansion and reference temperature must be finite.");
        }
        std::array<long double, 4> samples{};
        long double mean = 0;
        for (std::size_t a = 0; a < 4; ++a) {
            const auto value = temperature[impl_->cells[e].nodes[a]];
            if (!std::isfinite(value)) throw std::invalid_argument("Non-finite nodal temperature.");
            samples[a] = alpha * (static_cast<long double>(value) - reference);
            mean += samples[a] * 0.25L;
        }
        moments[e].mean = static_cast<double>(mean);
        for (const auto value : samples) {
            const long double difference = value - mean;
            moments[e].variance += difference * difference / 20.0L;
        }
    }
    return impl_->solve(moments, options, load_start);
}

LinearElasticity3DResult PreparedElasticity3D::solve_cell_eigenstrain(
    std::span<const double> strain, const linalg::SolverOptions& options
) {
    const auto load_start = Clock::now();
    if (!strain.empty() && strain.size() != impl_->cells.size()) {
        throw std::invalid_argument("Eigenstrain must cover every tetrahedron.");
    }
    std::vector<CellMoment> moments;
    moments.reserve(strain.size());
    for (const auto value : strain) {
        if (!std::isfinite(value)) throw std::invalid_argument("Non-finite cell eigenstrain.");
        moments.push_back({value, 0});
    }
    return impl_->solve(moments, options, load_start);
}

LinearElasticity3DResult solve_linear_elasticity_3d(
    const LinearElasticity3DProblem& problem, const linalg::SolverOptions& options
) {
    const auto& alpha = problem.cell_thermal_expansion_per_kelvin;
    const auto& delta = problem.cell_temperature_change_kelvin;
    if (alpha.size() != delta.size() || (!alpha.empty() && alpha.size() != problem.mesh.elements.size())) {
        throw std::invalid_argument("3D thermal strain requires matching per-cell expansion and temperature arrays.");
    }
    std::vector<double> strain(alpha.size());
    for (std::size_t e = 0; e < strain.size(); ++e) {
        if (!std::isfinite(alpha[e]) || !std::isfinite(delta[e])) {
            throw std::invalid_argument("Thermal expansion and temperature changes must be finite.");
        }
        strain[e] = alpha[e] * delta[e];
        if (!std::isfinite(strain[e])) throw std::overflow_error("Non-finite 3D thermal eigenstrain.");
    }
    PreparedElasticity3D prepared(problem.mesh, problem.cell_materials, problem.prescribed_displacements,
                                  problem.cell_body_force_density, problem.face_tractions,
                                  problem.component_constraints);
    return prepared.solve_cell_eigenstrain(strain, options);
}

} // namespace nova::physics
