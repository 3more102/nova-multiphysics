#include "nova/fem/conservative_face_heat_3d.hpp"

#include "nova/fem/energy_balance_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"
#include "nova/fem/tetrahedron_h1.hpp"
#include "nova/linalg/conjugate_gradient.hpp"

#include "p1_integration_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using FaceKey = BoundaryFace3D;
using Vector3L = std::array<long double, 3>;

struct FaceRecord {
    ReconciledFaceHeat3D face;
    std::size_t opposite{};
    Vector3L outward_area{};
    bool adjustable{};
};

[[nodiscard]] FaceKey canonical(FaceKey face) {
    std::sort(face.begin(), face.end());
    return face;
}

[[nodiscard]] Vector3L difference(const Point3D& a, const Point3D& b) {
    return {static_cast<long double>(a.x) - b.x,
            static_cast<long double>(a.y) - b.y,
            static_cast<long double>(a.z) - b.z};
}

[[nodiscard]] Vector3L cross(const Vector3L& a, const Vector3L& b) {
    return {a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]};
}

[[nodiscard]] long double dot(const Vector3L& a, const Vector3L& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

[[nodiscard]] long double area_length(const Vector3L& v) {
    return std::hypot(v[0], v[1], v[2]);
}

[[nodiscard]] double checked(long double x, const char* what) {
    if (!std::isfinite(x) ||
        std::abs(x) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error(what);
    }
    return static_cast<double>(x);
}

struct DisjointSet {
    explicit DisjointSet(std::size_t n) : parent(n) {
        std::iota(parent.begin(), parent.end(), std::size_t{0});
    }
    std::size_t find(std::size_t a) {
        if (parent[a] != a) parent[a] = find(parent[a]);
        return parent[a];
    }
    void join(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a != b) parent[b] = a;
    }
    std::vector<std::size_t> parent;
};

// The same degree-two triangle quadrature used by the Galerkin solve.
[[nodiscard]] long double imposed_neumann(
    const TetrahedralMesh3D& mesh, const NeumannBoundaryFace3D& bc
) {
    const double jacobian = detail::face_jacobian(mesh, bc.face);
    long double rate = 0.0L;
    for (const auto& q : triangle_face_degree_two_quadrature) {
        const auto p = detail::simplex_point(mesh, bc.face, q.shape);
        rate += static_cast<long double>(
            bc.outward_flux(p.x, p.y, p.z)
        ) * jacobian * q.weight;
    }
    return rate;
}

[[nodiscard]] long double imposed_robin(
    const TetrahedralMesh3D& mesh, const RobinBoundaryFace3D& bc,
    std::span<const double> values
) {
    const double jacobian = detail::face_jacobian(mesh, bc.face);
    long double rate = 0.0L;
    for (const auto& q : triangle_face_degree_two_quadrature) {
        const auto p = detail::simplex_point(mesh, bc.face, q.shape);
        const double h = bc.transfer_coefficient(p.x, p.y, p.z);
        const double ambient = bc.ambient_value(p.x, p.y, p.z);
        long double temperature = 0.0L;
        for (std::size_t j = 0; j < 3; ++j) {
            temperature += static_cast<long double>(values[bc.face[j]]) * q.shape[j];
        }
        rate += static_cast<long double>(h) *
                (temperature - ambient) * jacobian * q.weight;
    }
    return rate;
}

[[nodiscard]] std::vector<long double> element_generation(
    const ScalarDiffusion3DProblem& problem
) {
    const auto& mesh = problem.mesh;
    std::vector<long double> generation(mesh.elements.size(), 0.0L);
    for (std::size_t cell = 0; cell < mesh.elements.size(); ++cell) {
        const auto& tetra = mesh.elements[cell];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, cell);
        const long double jacobian =
            std::abs(static_cast<long double>(basis.signed_jacobian_determinant));
        const double cell_source = problem.cell_source.empty()
            ? 0.0 : problem.cell_source[cell];
        for (const auto& q : tetrahedron_degree_two_quadrature) {
            const auto p = detail::simplex_point(mesh, tetra, q.shape);
            const double spatial = problem.source(p.x, p.y, p.z);
            const double density = spatial + cell_source;
            if (!std::isfinite(density)) {
                throw std::overflow_error("Reconstruction source density overflowed.");
            }
            generation[cell] += static_cast<long double>(density) *
                                jacobian * q.weight;
        }
        (void)checked(generation[cell], "Reconstruction heat generation overflowed.");
    }
    return generation;
}

} // namespace

FaceHeatReconciliation3D reconcile_p1_face_heat_rates_3d(
    const ScalarDiffusion3DProblem& problem,
    const std::span<const double> nodal_values,
    const linalg::SolverOptions& correction_solver_options
) {
    if (!problem.dirichlet_constraints.empty()) {
        throw std::invalid_argument(
            "Conservative face rates require face Dirichlet data; "
            "pointwise Dirichlet reactions cannot be distributed to faces."
        );
    }

    // Validates all face assignments, constrained temperatures, source data,
    // quadrature callbacks and material representation using the steady solver's
    // existing independent energy audit. No solver is invoked.
    (void)evaluate_diffusion_energy_balance_3d(problem, nodal_values);
    const auto raw_flux = evaluate_p1_tetrahedral_flux_3d(problem, nodal_values);
    const auto& mesh = problem.mesh;
    const std::size_t n = mesh.elements.size();

    std::map<FaceKey, FaceRecord> faces;
    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto& tetra = mesh.elements[cell];
        for (std::size_t missing = 0; missing < 4; ++missing) {
            FaceKey key{};
            std::size_t next = 0;
            for (std::size_t local = 0; local < 4; ++local) {
                if (local != missing) key[next++] = tetra[local];
            }
            key = canonical(key);
            const auto [it, fresh] = faces.try_emplace(key);
            if (fresh) {
                it->second.face.nodes = key;
                it->second.face.first_cell = cell;
                it->second.opposite = tetra[missing];
            } else {
                it->second.face.interior = true;
                it->second.face.second_cell = cell;
            }
        }
    }

    std::set<FaceKey> dirichlet;
    std::map<FaceKey, long double> neumann;
    std::map<FaceKey, long double> robin;
    for (const auto& bc : problem.dirichlet_faces) {
        dirichlet.insert(canonical(bc.face));
    }
    for (const auto& bc : problem.neumann_faces) {
        neumann.emplace(canonical(bc.face), imposed_neumann(mesh, bc));
    }
    for (const auto& bc : problem.robin_faces) {
        robin.emplace(canonical(bc.face), imposed_robin(mesh, bc, nodal_values));
    }

    DisjointSet components(n);
    for (const auto& [key, record] : faces) {
        (void)key;
        if (record.face.interior) {
            components.join(record.face.first_cell, record.face.second_cell);
        }
    }

    std::vector<bool> anchored(n, false);
    std::vector<long double> outgoing_baseline(n, 0.0L);
    // Compatibility uses only generation and signed prescribed exterior rates;
    // internal face fluxes cancel and must never set the component scale.
    std::vector<long double> exterior_baseline_rate(n, 0.0L);
    std::vector<long double> exterior_baseline_magnitude(n, 0.0L);
    for (auto& [key, record] : faces) {
        const auto& a = mesh.nodes[key[0]];
        const auto ab = difference(mesh.nodes[key[1]], a);
        const auto ac = difference(mesh.nodes[key[2]], a);
        const auto to_opposite = difference(mesh.nodes[record.opposite], a);
        auto area_vector = cross(ab, ac);
        for (auto& coordinate : area_vector) coordinate *= 0.5L;
        if (dot(area_vector, to_opposite) > 0.0L) {
            for (auto& coordinate : area_vector) coordinate = -coordinate;
        }
        const long double area = area_length(area_vector);
        record.outward_area = area_vector;
        record.face.area = checked(area, "Face area overflows.");
        if (!(record.face.area > 0.0)) {
            throw std::invalid_argument("Reconstruction requires nondegenerate faces.");
        }

        const auto& first = raw_flux[record.face.first_cell].flux;
        const Vector3L first_vector{first.x, first.y, first.z};
        long double baseline = dot(first_vector, area_vector);
        if (record.face.interior) {
            const auto& second = raw_flux[record.face.second_cell].flux;
            const Vector3L second_vector{second.x, second.y, second.z};
            baseline = 0.5L * (baseline + dot(second_vector, area_vector));
            record.face.condition = FaceHeatCondition3D::Interior;
            record.adjustable = true;
        } else if (dirichlet.contains(key)) {
            record.face.condition = FaceHeatCondition3D::Dirichlet;
            record.adjustable = true;
            anchored[components.find(record.face.first_cell)] = true;
        } else if (const auto it = neumann.find(key); it != neumann.end()) {
            record.face.condition = FaceHeatCondition3D::Neumann;
            record.adjustable = false;
            record.face.outward_heat = checked(
                it->second, "Prescribed Neumann rate overflows."
            );
        } else if (const auto it = robin.find(key); it != robin.end()) {
            record.face.condition = FaceHeatCondition3D::Robin;
            record.adjustable = false;
            record.face.outward_heat = checked(
                it->second, "Prescribed Robin rate overflows."
            );
        } else {
            record.face.condition = FaceHeatCondition3D::Insulated;
            record.adjustable = false;
            record.face.outward_heat = 0.0;
        }
        record.face.raw_outward_heat =
            checked(baseline, "Raw integrated face heat rate overflows.");
        if (record.adjustable) {
            record.face.outward_heat = record.face.raw_outward_heat;
        }
        const long double heat = record.face.outward_heat;
        outgoing_baseline[record.face.first_cell] += heat;
        if (record.face.interior) {
            outgoing_baseline[record.face.second_cell] -= heat;
        } else {
            exterior_baseline_rate[record.face.first_cell] += heat;
            exterior_baseline_magnitude[record.face.first_cell] += std::abs(heat);
        }
    }

    const auto generation = element_generation(problem);
    std::vector<long double> imbalance(n);
    for (std::size_t i = 0; i < n; ++i) {
        imbalance[i] = generation[i] - outgoing_baseline[i];
        (void)checked(imbalance[i], "Reconstruction cell imbalance overflows.");
    }

    // The incidence graph's Laplacian has one null mode per connected
    // component without a Dirichlet face. Check compatibility *before*
    // grounding one multiplier in each such component.
    std::vector<long double> component_defect(n, 0.0L);
    std::vector<long double> component_rate_scale(n, 0.0L);
    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto root = components.find(cell);
        component_defect[root] += generation[cell] - exterior_baseline_rate[cell];
        component_rate_scale[root] += std::abs(generation[cell]) +
                                      exterior_baseline_magnitude[cell];
    }

    std::vector<bool> gauges(n, false);
    std::vector<bool> seen(n, false);
    std::size_t floating_count = 0;
    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto root = components.find(cell);
        if (anchored[root] || seen[root]) continue;
        seen[root] = true;
        const long double limit = 1.0e-10L * component_rate_scale[root];
        if (std::abs(component_defect[root]) > limit) {
            throw std::invalid_argument(
                "Unanchored component has incompatible net heat generation "
                "and prescribed exterior heat flow."
            );
        }
        gauges[cell] = true;
        ++floating_count;
    }

    // Scale each disconnected component independently. Without this, one
    // large-power component can set CG's global residual scale and suppress
    // corrections needed by a physically valid, much smaller component.
    std::vector<long double> component_rhs_scale(n, 0.0L);
    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto root = components.find(cell);
        component_rhs_scale[root] = std::max(
            component_rhs_scale[root], std::abs(imbalance[cell])
        );
    }
    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto root = components.find(cell);
        if (component_rhs_scale[root] == 0.0L) {
            component_rhs_scale[root] = 1.0L;
        }
    }

    std::vector<linalg::Triplet> entries;
    if (n > entries.max_size() || faces.size() > (entries.max_size() - n) / 4U) {
        throw std::length_error("Face reconciliation matrix exceeds storage limits.");
    }
    entries.reserve(faces.size() * 4U + n);
    for (const auto& [key, record] : faces) {
        (void)key;
        if (!record.adjustable) continue;
        // Minimize the physical surface L2 correction:
        // integral_F (delta_heat_rate / area_F)^2 dA = delta_heat_rate^2 / area_F.
        const double weight = record.face.area;
        if (!std::isfinite(weight) || weight <= 0.0) {
            throw std::overflow_error("Face reconciliation weight is invalid.");
        }
        const auto a = record.face.first_cell;
        if (!gauges[a]) entries.push_back({a, a, weight});
        if (record.face.interior) {
            const auto b = record.face.second_cell;
            if (!gauges[b]) entries.push_back({b, b, weight});
            if (!gauges[a] && !gauges[b]) {
                entries.push_back({a, b, -weight});
                entries.push_back({b, a, -weight});
            }
        }
    }
    std::vector<double> rhs(n, 0.0);
    for (std::size_t cell = 0; cell < n; ++cell) {
        if (gauges[cell]) {
            entries.push_back({cell, cell, 1.0});
        } else {
            rhs[cell] = checked(
                imbalance[cell] / component_rhs_scale[components.find(cell)],
                "Scaled cell correction RHS overflows."
            );
        }
    }
    const auto matrix = linalg::CsrMatrix::from_triplets(n, n, entries);
    auto multipliers = linalg::solve_conjugate_gradient(
        matrix, rhs, correction_solver_options
    );
    if (!multipliers.converged) {
        throw std::runtime_error(
            "Face heat reconciliation CG did not converge: " + multipliers.message
        );
    }

    FaceHeatReconciliation3D result;
    result.multiplier_solver = std::move(multipliers);
    result.floating_components = floating_count;
    result.generated_heat_per_cell.reserve(n);
    result.cell_balance_defect.resize(n, 0.0);
    std::vector<long double> outgoing(n, 0.0L);
    std::vector<long double> outgoing_magnitude(n, 0.0L);
    long double correction_density_square = 0.0L;
    long double maximum_correction = 0.0L;
    long double total_exterior = 0.0L;
    long double total_generated = 0.0L;

    for (auto& [key, record] : faces) {
        (void)key;
        long double correction_magnitude = 0.0L;
        if (record.adjustable) {
            const auto a = record.face.first_cell;
            long double multiplier_difference = result.multiplier_solver.solution[a];
            if (record.face.interior) {
                multiplier_difference -=
                    result.multiplier_solver.solution[record.face.second_cell];
            }
            const long double correction =
                static_cast<long double>(record.face.area) * multiplier_difference *
                component_rhs_scale[components.find(a)];
            correction_magnitude = std::abs(correction);
            record.face.outward_heat = checked(
                static_cast<long double>(record.face.raw_outward_heat) +
                correction, "Reconciled face rate overflows."
            );
            maximum_correction = std::max(maximum_correction, std::abs(correction));
            const long double density = correction / record.face.area;
            correction_density_square +=
                static_cast<long double>(record.face.area) * density * density;
        }
        const long double heat = record.face.outward_heat;
        // Include the magnitudes of both terms in an adjustable-face
        // subtraction. If a large raw flux is cancelled by a large correction,
        // the remaining roundoff must be judged against those operation scales.
        const long double face_magnitude = record.adjustable
            ? std::abs(static_cast<long double>(record.face.raw_outward_heat)) +
                  correction_magnitude
            : std::abs(heat);
        outgoing[record.face.first_cell] += heat;
        outgoing_magnitude[record.face.first_cell] += face_magnitude;
        if (record.face.interior) {
            outgoing[record.face.second_cell] -= heat;
            outgoing_magnitude[record.face.second_cell] += face_magnitude;
        } else {
            total_exterior += heat;
        }
        result.faces.push_back(record.face);
    }

    for (std::size_t cell = 0; cell < n; ++cell) {
        const long double defect = generation[cell] - outgoing[cell];
        // Use the sum of absolute face rates, not the net outgoing sum:
        // large Neumann/Robin rates may cancel while the local defect is tiny.
        const long double scale =
            std::abs(generation[cell]) + outgoing_magnitude[cell];
        if (std::abs(defect) > 1.0e-8L * scale) {
            throw std::runtime_error(
                "Face heat correction has an excessive local conservation defect."
            );
        }
        result.generated_heat_per_cell.push_back(
            checked(generation[cell], "Cell generation exceeds finite range.")
        );
        result.cell_balance_defect[cell] = checked(
            defect, "Cell conservation defect overflows."
        );
        result.max_absolute_cell_defect = std::max(
            result.max_absolute_cell_defect, std::abs(result.cell_balance_defect[cell])
        );
        total_generated += generation[cell];
    }
    result.max_absolute_correction =
        checked(maximum_correction, "Face flux correction overflows.");
    result.correction_l2_flux_density = checked(
        std::sqrt(correction_density_square), "Face correction norm overflows."
    );
    result.total_generated_heat = checked(total_generated, "Total generation overflows.");
    result.total_exterior_outward_heat = checked(
        total_exterior, "Total exterior heat overflows."
    );
    return result;
}

} // namespace nova::fem
