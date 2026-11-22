#include "nova/physics/linear_elasticity_3d.hpp"

#include "nova/fem/tetrahedron_h1.hpp"
#include "constrained_diffusion.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::physics {
namespace {

struct LameParameters {
    double lambda{};
    double mu{};
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
    if (!std::isfinite(mu) || !(mu > 0.0) || !std::isfinite(lambda) ||
        !std::isfinite(lambda + 2.0 * mu / 3.0) ||
        !(lambda + 2.0 * mu / 3.0 > 0.0)) {
        throw std::overflow_error("3D elasticity Lame coefficients are not representable.");
    }
    return {lambda, mu};
}

[[nodiscard]] bool finite(const ElasticVector3D& x) {
    return std::all_of(x.begin(), x.end(),
                       [](double value) { return std::isfinite(value); });
}

[[nodiscard]] std::size_t find_root(
    std::vector<std::size_t>& parents, std::size_t node
) {
    while (parents[node] != node) {
        parents[node] = parents[parents[node]];
        node = parents[node];
    }
    return node;
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

// Full-vector vertex constraints remove all six rigid motions if, and only
// if, the anchored points contain a noncollinear triple on each component.
void require_rigid_mode_anchors(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const DisplacementConstraint3D> constraints
) {
    std::vector<std::size_t> parents(mesh.nodes.size());
    for (std::size_t i = 0; i < parents.size(); ++i) parents[i] = i;
    for (const auto& cell : mesh.elements) {
        const auto first = find_root(parents, cell[0]);
        for (std::size_t i = 1; i < cell.size(); ++i) {
            const auto other = find_root(parents, cell[i]);
            if (first != other) parents[other] = first;
        }
    }
    std::map<std::size_t, std::vector<std::size_t>> fixed_by_component;
    for (const auto& constraint : constraints) {
        fixed_by_component[find_root(parents, constraint.node)].push_back(
            constraint.node);
    }
    std::set<std::size_t> components;
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
        components.insert(find_root(parents, i));
    }
    for (const auto component : components) {
        const auto found = fixed_by_component.find(component);
        if (found == fixed_by_component.end() || found->second.size() < 3) {
            throw std::invalid_argument(
                "Each 3D elastic component needs three noncollinear fully fixed vertices.");
        }
        const auto& nodes = found->second;
        const auto& first = mesh.nodes[nodes.front()];
        ElasticVector3D direction{};
        double longest = 0.0;
        for (const auto node : nodes) {
            const auto trial = between(first, mesh.nodes[node]);
            const double length2 = norm_squared(trial);
            if (length2 > longest) {
                longest = length2;
                direction = trial;
            }
        }
        if (!(longest > 0.0) || !std::isfinite(longest)) {
            throw std::invalid_argument("Degenerate 3D elastic displacement anchors.");
        }
        bool independent = false;
        for (const auto node : nodes) {
            const auto trial = between(first, mesh.nodes[node]);
            const double trial_norm2 = norm_squared(trial);
            if (trial_norm2 == 0.0) continue;
            const double sin2_numerator = norm_squared(cross(direction, trial));
            if (std::isfinite(sin2_numerator) &&
                sin2_numerator > 1.0e-24 * longest * trial_norm2) {
                independent = true;
                break;
            }
        }
        if (!independent) {
            throw std::invalid_argument(
                "Collinear 3D elastic displacement anchors leave a rigid rotation.");
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
    const double volumetric = mat.lambda * (eps[0] + eps[1] + eps[2]);
    return {
        volumetric + 2.0 * mat.mu * eps[0],
        volumetric + 2.0 * mat.mu * eps[1],
        volumetric + 2.0 * mat.mu * eps[2],
        2.0 * mat.mu * eps[3],
        2.0 * mat.mu * eps[4],
        2.0 * mat.mu * eps[5]
    };
}

} // namespace

LinearElasticity3DResult solve_linear_elasticity_3d(
    const LinearElasticity3DProblem& problem,
    const linalg::SolverOptions& solver_options
) {
    problem.mesh.validate();
    const auto node_count = problem.mesh.nodes.size();
    const auto cell_count = problem.mesh.elements.size();
    if (node_count > std::numeric_limits<std::size_t>::max() / 3U) {
        throw std::length_error("3D elasticity DOF count overflows.");
    }
    if (problem.cell_materials.size() != cell_count ||
        (!problem.cell_body_force_density.empty() &&
         problem.cell_body_force_density.size() != cell_count)) {
        throw std::invalid_argument(
            "3D elasticity requires one material and optional body force per cell.");
    }
    std::vector<LameParameters> materials;
    materials.reserve(cell_count);
    for (const auto& material : problem.cell_materials) {
        materials.push_back(lame(material));
    }
    for (const auto& f : problem.cell_body_force_density) {
        if (!finite(f)) {
            throw std::invalid_argument("3D elastic body-force density must be finite.");
        }
    }

    std::vector<std::optional<double>> prescribed(3U * node_count);
    std::set<std::size_t> fixed_vertices;
    for (const auto& condition : problem.prescribed_displacements) {
        if (condition.node >= node_count) {
            throw std::out_of_range("3D elastic constraint references a missing vertex.");
        }
        if (!fixed_vertices.insert(condition.node).second || !finite(condition.value)) {
            throw std::invalid_argument(
                "Duplicate 3D elastic vertex constraint or non-finite displacement.");
        }
        for (std::size_t i = 0; i < 3U; ++i) {
            prescribed[3U * condition.node + i] = condition.value[i];
        }
    }
    require_rigid_mode_anchors(problem.mesh, problem.prescribed_displacements);

    std::set<fem::BoundaryFace3D> exterior;
    for (const auto& face : extract_boundary_faces_3d(problem.mesh)) {
        exterior.insert(sorted_face(face));
    }
    std::set<fem::BoundaryFace3D> assigned_tractions;
    for (const auto& condition : problem.face_tractions) {
        if (!exterior.contains(sorted_face(condition.face)) ||
            !assigned_tractions.insert(sorted_face(condition.face)).second ||
            !finite(condition.traction)) {
            throw std::invalid_argument(
                "3D elastic traction must be finite and assigned once to an exterior face.");
        }
    }

    constexpr std::size_t entries_per_cell = 12U * 12U;
    std::vector<linalg::Triplet> entries;
    if (cell_count > entries.max_size() / entries_per_cell) {
        throw std::length_error("3D elasticity stiffness exceeds storage capacity.");
    }
    entries.reserve(cell_count * entries_per_cell);
    std::vector<double> rhs(3U * node_count, 0.0);

    for (std::size_t e = 0; e < cell_count; ++e) {
        const auto& tet = problem.mesh.elements[e];
        const auto geometry = fem::evaluate_tetrahedron_p1_physical_gradients_3d(
            problem.mesh, e);
        const double volume = std::abs(geometry.signed_jacobian_determinant) / 6.0;
        const auto& mat = materials[e];
        const auto& grad = geometry.gradients;
        const auto force = problem.cell_body_force_density.empty()
            ? ElasticVector3D{} : problem.cell_body_force_density[e];
        for (std::size_t a = 0; a < 4U; ++a) {
            for (std::size_t i = 0; i < 3U; ++i) {
                const auto row = 3U * tet[a] + i;
                rhs[row] += volume * force[i] / 4.0;
                for (std::size_t b = 0; b < 4U; ++b) {
                    const double dot = grad[a][0] * grad[b][0] +
                                       grad[a][1] * grad[b][1] +
                                       grad[a][2] * grad[b][2];
                    for (std::size_t j = 0; j < 3U; ++j) {
                        const double stiffness = volume * (
                            mat.lambda * grad[a][i] * grad[b][j] +
                            mat.mu * grad[a][j] * grad[b][i] +
                            (i == j ? mat.mu * dot : 0.0));
                        if (!std::isfinite(stiffness)) {
                            throw std::overflow_error("Non-finite 3D elastic stiffness.");
                        }
                        entries.push_back({row, 3U * tet[b] + j, stiffness});
                    }
                }
                if (!std::isfinite(rhs[row])) {
                    throw std::overflow_error("Non-finite 3D elastic nodal load.");
                }
            }
        }
    }
    for (const auto& condition : problem.face_tractions) {
        const auto& face = condition.face;
        const auto& p0 = problem.mesh.nodes[face[0]];
        const auto edge1 = between(p0, problem.mesh.nodes[face[1]]);
        const auto edge2 = between(p0, problem.mesh.nodes[face[2]]);
        const double area = 0.5 * std::sqrt(norm_squared(cross(edge1, edge2)));
        if (!(area > 0.0) || !std::isfinite(area)) {
            throw std::overflow_error("Invalid 3D elastic traction face area.");
        }
        for (const auto node : face) {
            for (std::size_t i = 0; i < 3U; ++i) {
                const auto row = 3U * node + i;
                rhs[row] += area * condition.traction[i] / 3.0;
                if (!std::isfinite(rhs[row])) {
                    throw std::overflow_error("Non-finite 3D elastic traction load.");
                }
            }
        }
    }

    auto solved = fem::detail::solve_constrained_diffusion(
        std::move(entries), std::move(rhs), prescribed, solver_options);

    LinearElasticity3DResult result;
    result.linear_solver = std::move(solved);
    result.displacement.resize(node_count);
    for (std::size_t node = 0; node < node_count; ++node) {
        for (std::size_t i = 0; i < 3U; ++i) {
            result.displacement[node][i] = result.linear_solver.solution[3U * node + i];
        }
    }
    result.cell_strain.reserve(cell_count);
    result.cell_stress.reserve(cell_count);
    for (std::size_t e = 0; e < cell_count; ++e) {
        const auto geometry = fem::evaluate_tetrahedron_p1_physical_gradients_3d(
            problem.mesh, e);
        const double volume = std::abs(geometry.signed_jacobian_determinant) / 6.0;
        std::array<ElasticVector3D, 3> gradient{};
        for (std::size_t a = 0; a < 4U; ++a) {
            const auto node = problem.mesh.elements[e][a];
            for (std::size_t i = 0; i < 3U; ++i) {
                for (std::size_t j = 0; j < 3U; ++j) {
                    gradient[i][j] += result.displacement[node][i] *
                                      geometry.gradients[a][j];
                }
            }
        }
        const auto strain = strain_from_gradient(gradient);
        const auto stress = stress_from_strain(strain, materials[e]);
        double work = 0.0;
        for (std::size_t i = 0; i < 6U; ++i) {
            work += (i < 3U ? 1.0 : 2.0) * strain[i] * stress[i];
        }
        result.strain_energy_joules += 0.5 * volume * work;
        if (!std::isfinite(result.strain_energy_joules)) {
            throw std::overflow_error("Non-finite 3D elastic strain energy.");
        }
        result.cell_strain.push_back(strain);
        result.cell_stress.push_back(stress);
    }
    return result;
}

} // namespace nova::physics
