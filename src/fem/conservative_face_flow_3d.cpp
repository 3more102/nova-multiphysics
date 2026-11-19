#include "nova/fem/conservative_face_flow_3d.hpp"

#include "nova/fem/energy_balance_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"
#include "nova/fem/simplex_quadrature.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

#include "p1_integration_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using FaceKey = BoundaryFace3D;
constexpr std::size_t missing = std::numeric_limits<std::size_t>::max();

struct WorkingFace {
    ConservativeFaceFlow3D report;
    Point3D area_vector{}; // Outward face area vector of owner_cell.
    long double rate{};
    long double raw{};
};

FaceKey key(FaceKey vertices) {
    std::sort(vertices.begin(), vertices.end());
    return vertices;
}

Point3D minus(const Point3D& a, const Point3D& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Point3D cross(const Point3D& a, const Point3D& b) {
    return {
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

long double dot(const Point3D& a, const Point3D& b) {
    return static_cast<long double>(a.x)*b.x +
           static_cast<long double>(a.y)*b.y +
           static_cast<long double>(a.z)*b.z;
}

Point3D outward_area(const TetrahedralMesh3D& mesh,
                     const FaceKey& face, std::size_t opposite) {
    const auto& a = mesh.nodes[face[0]];
    const auto area = cross(minus(mesh.nodes[face[1]], a),
                            minus(mesh.nodes[face[2]], a));
    const double orientation = static_cast<double>(dot(area, minus(mesh.nodes[opposite], a)));
    const double sign = orientation > 0.0 ? -0.5 : 0.5;
    const Point3D result{sign*area.x, sign*area.y, sign*area.z};
    if (!std::isfinite(result.x) || !std::isfinite(result.y) || !std::isfinite(result.z) ||
        orientation == 0.0) {
        throw std::invalid_argument("Cannot orient a conservative tetrahedral face.");
    }
    return result;
}

double checked(long double v, const char* message) {
    if (!std::isfinite(v) ||
        std::abs(v) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error(message);
    }
    return static_cast<double>(v);
}

template<class Integrand>
long double integrate_face(const TetrahedralMesh3D& mesh,
                           const FaceKey& face, Integrand&& value) {
    const double jac = detail::face_jacobian(mesh, face);
    long double sum = 0.0L;
    for (const auto& q : triangle_face_degree_two_quadrature) {
        const Point3D p = detail::simplex_point(mesh, face, q.shape);
        const double sample = value(p, q.shape);
        if (!std::isfinite(sample)) {
            throw std::invalid_argument("Non-finite boundary flow in conservation projection.");
        }
        sum += static_cast<long double>(sample) * jac * q.weight;
    }
    return sum;
}

} // namespace

ConservativeHeatFlow3D reconstruct_conservative_heat_flow_3d(
    const ScalarDiffusion3DProblem& problem,
    std::span<const double> nodal_values
) {
    // This also checks finite data, physical material, mutually exclusive
    // boundary assignments, prescribed temperatures, and source callbacks.
    (void)evaluate_diffusion_energy_balance_3d(problem, nodal_values);
    const auto& mesh = problem.mesh;
    const auto fields = evaluate_p1_tetrahedral_flux_3d(problem, nodal_values);
    const std::size_t n = mesh.elements.size();

    std::map<FaceKey, std::size_t> lookup;
    std::vector<WorkingFace> faces;
    faces.reserve(n * 2U + 2U);
    constexpr std::array<std::array<std::size_t, 3>, 4> local_faces{{
        {{1,2,3}}, {{0,2,3}}, {{0,1,3}}, {{0,1,2}}
    }};
    for (std::size_t c = 0; c < n; ++c) {
        const auto& tetra = mesh.elements[c];
        for (std::size_t opposite = 0; opposite < 4; ++opposite) {
            const auto& ids = local_faces[opposite];
            const FaceKey vertices{tetra[ids[0]], tetra[ids[1]], tetra[ids[2]]};
            const auto face_key = key(vertices);
            const auto found = lookup.find(face_key);
            if (found == lookup.end()) {
                const auto area = outward_area(mesh, vertices, tetra[opposite]);
                const long double raw = dot(fields[c].flux, area);
                const std::size_t index = faces.size();
                lookup.emplace(face_key, index);
                WorkingFace entry;
                entry.report.nodes = face_key;
                entry.report.owner_cell = c;
                entry.area_vector = area;
                entry.rate = raw;
                entry.raw = raw;
                faces.push_back(entry);
            } else {
                auto& face = faces[found->second];
                if (face.report.neighbor_cell) {
                    throw std::invalid_argument("Non-manifold face in conservative heat flow.");
                }
                face.report.neighbor_cell = c;
                // The shared rate is continuous by construction. The normal
                // on both sides is the owner's outward area vector.
                face.rate = 0.5L * (face.raw + dot(fields[c].flux, face.area_vector));
                face.raw = face.rate;
            }
        }
    }

    auto face_for = [&](FaceKey vertices) -> WorkingFace& {
        const auto found = lookup.find(key(vertices));
        if (found == lookup.end() || faces[found->second].report.neighbor_cell) {
            throw std::invalid_argument("Boundary condition requires an exterior face.");
        }
        return faces[found->second];
    };
    for (const auto& condition : problem.dirichlet_faces) {
        auto& face = face_for(condition.face);
        face.report.kind = FaceFlowBoundaryKind3D::dirichlet;
    }
    for (const auto& condition : problem.neumann_faces) {
        auto& face = face_for(condition.face);
        face.report.kind = FaceFlowBoundaryKind3D::neumann;
        face.rate = integrate_face(mesh, condition.face, [&](const Point3D& p, const auto&) {
            return condition.outward_flux(p.x, p.y, p.z);
        });
        face.raw = face.rate;
    }
    for (const auto& condition : problem.robin_faces) {
        auto& face = face_for(condition.face);
        face.report.kind = FaceFlowBoundaryKind3D::robin;
        face.rate = integrate_face(mesh, condition.face, [&](const Point3D& p, const auto& shape) {
            const double h = condition.transfer_coefficient(p.x,p.y,p.z);
            const double ambient = condition.ambient_value(p.x,p.y,p.z);
            long double temp = 0.0L;
            for (std::size_t k = 0; k < 3; ++k) {
                temp += static_cast<long double>(shape[k]) * nodal_values[condition.face[k]];
            }
            const long double flow = static_cast<long double>(h) * (temp - ambient);
            return checked(flow, "Robin boundary flow exceeds finite range.");
        });
        face.raw = face.rate;
    }
    for (auto& face : faces) {
        if (!face.report.neighbor_cell &&
            face.report.kind == FaceFlowBoundaryKind3D::interior) {
            face.report.kind = FaceFlowBoundaryKind3D::insulated;
            face.rate = 0.0L;
            face.raw = 0.0L;
        }
    }

    std::vector<long double> generation(n, 0.0L);
    for (std::size_t c = 0; c < n; ++c) {
        const auto& tetra = mesh.elements[c];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, c);
        const double jac = std::abs(basis.signed_jacobian_determinant);
        const double cell_source = problem.cell_source.empty() ? 0.0 : problem.cell_source[c];
        for (const auto& q : tetrahedron_degree_two_quadrature) {
            const auto p = detail::simplex_point(mesh, tetra, q.shape);
            const double spatial = problem.source(p.x,p.y,p.z);
            const double density = spatial + cell_source;
            if (!std::isfinite(density)) {
                throw std::overflow_error("Combined heat source exceeds finite range.");
            }
            generation[c] += static_cast<long double>(density)*jac*q.weight;
        }
    }

    std::vector<std::vector<std::pair<std::size_t,std::size_t>>> adjacent(n);
    std::vector<std::size_t> outlet(n, missing);
    auto defects = generation;
    for (std::size_t f = 0; f < faces.size(); ++f) {
        const auto& face = faces[f];
        const auto owner = face.report.owner_cell;
        defects[owner] -= face.rate;
        if (face.report.neighbor_cell) {
            const auto other = *face.report.neighbor_cell;
            defects[other] += face.rate;
            adjacent[owner].push_back({other, f});
            adjacent[other].push_back({owner, f});
        } else if (face.report.kind == FaceFlowBoundaryKind3D::dirichlet) {
            if (outlet[owner] == missing) outlet[owner] = f;
        }
    }
    const auto before = defects;

    std::vector<bool> component_visited(n, false);
    std::vector<std::size_t> parent(n, missing);
    std::vector<std::size_t> via(n, missing);
    for (std::size_t start = 0; start < n; ++start) {
        if (component_visited[start]) continue;
        std::vector<std::size_t> component{start};
        component_visited[start] = true;
        for (std::size_t k = 0; k < component.size(); ++k) {
            for (const auto& [neighbor, edge] : adjacent[component[k]]) {
                (void)edge;
                if (!component_visited[neighbor]) {
                    component_visited[neighbor] = true;
                    component.push_back(neighbor);
                }
            }
        }
        std::size_t root = component.front();
        for (auto c : component) {
            if (outlet[c] != missing) { root = c; break; }
        }
        std::vector<std::size_t> order{root};
        parent[root] = root;
        for (std::size_t k = 0; k < order.size(); ++k) {
            const auto c = order[k];
            for (const auto& [neighbor, edge] : adjacent[c]) {
                if (parent[neighbor] == missing) {
                    parent[neighbor] = c;
                    via[neighbor] = edge;
                    order.push_back(neighbor);
                }
            }
        }
        for (std::size_t k = order.size(); k-- > 1;) {
            const std::size_t child = order[k];
            const std::size_t edge = via[child];
            const long double adjustment = defects[child];
            faces[edge].rate += faces[edge].report.owner_cell == child
                ? adjustment : -adjustment;
            defects[parent[child]] += adjustment;
            defects[child] = 0.0L;
        }
        if (outlet[root] != missing) {
            faces[outlet[root]].rate += defects[root];
            defects[root] = 0.0L;
        } else {
            long double scale = 1.0L;
            for (const auto c : component) scale += std::abs(generation[c]);
            for (const auto& f : faces) {
                if (!f.report.neighbor_cell &&
                    parent[f.report.owner_cell] != missing) {
                    scale += std::abs(f.rate);
                }
            }
            if (std::abs(defects[root]) > 1e-10L * scale) {
                throw std::invalid_argument(
                    "Incompatible source and prescribed boundary heat flows: "
                    "no exterior Dirichlet outlet on a connected cell component.");
            }
        }
    }

    ConservativeHeatFlow3D result;
    result.faces.reserve(faces.size());
    for (const auto& face : faces) {
        auto entry = face.report;
        entry.raw_outward_heat_rate = checked(face.raw, "Raw face heat rate overflow.");
        entry.outward_heat_rate = checked(face.rate, "Reconstructed face heat rate overflow.");
        result.maximum_absolute_face_correction =
            std::max(result.maximum_absolute_face_correction,
                     checked(std::abs(face.rate - face.raw), "Face correction overflow."));
        if (!entry.neighbor_cell) {
            result.total_exterior_outward_heat += entry.outward_heat_rate;
        }
        result.faces.push_back(entry);
    }
    // Recompute conservation defects from the exported double-precision rates,
    // not the long-double forest bookkeeping, so diagnostics describe the
    // actual face flows that downstream SDK consumers receive.
    std::vector<long double> exported_defects = generation;
    for (const auto& face : result.faces) {
        exported_defects[face.owner_cell] -= face.outward_heat_rate;
        if (face.neighbor_cell)
            exported_defects[*face.neighbor_cell] += face.outward_heat_rate;
    }
    result.generated_heat_per_cell.reserve(n);
    result.defect_before.reserve(n);
    result.defect_after.reserve(n);
    for (std::size_t c = 0; c < n; ++c) {
        const double g = checked(generation[c], "Cell generation overflow.");
        const double b = checked(before[c], "Unreconstructed cell balance overflow.");
        const double a = checked(exported_defects[c], "Reconstructed cell balance overflow.");
        result.generated_heat_per_cell.push_back(g);
        result.defect_before.push_back(b);
        result.defect_after.push_back(a);
        result.total_generated_heat += g;
        result.maximum_defect_before = std::max(result.maximum_defect_before, std::abs(b));
        result.maximum_defect_after = std::max(result.maximum_defect_after, std::abs(a));
    }
    return result;
}

} // namespace nova::fem
