#include "nova/fem/conservative_face_flux_3d.hpp"

#include "nova/fem/simplex_quadrature.hpp"
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

using Key = BoundaryFace3D;
constexpr std::size_t no_cell = std::numeric_limits<std::size_t>::max();

Key sorted(Key face) {
    std::sort(face.begin(), face.end());
    return face;
}

double checked(long double value, const char* message) {
    if (!std::isfinite(value) ||
        std::abs(value) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error(message);
    }
    return static_cast<double>(value);
}

struct Face {
    std::size_t cell{};
    std::size_t opposite{};
    std::size_t neighbor{no_cell};
    double area{};
    double rate{};
    bool adjustable{false};
};

class DisjointSet {
public:
    explicit DisjointSet(std::size_t n) : parent_(n) {
        std::iota(parent_.begin(), parent_.end(), std::size_t{0});
    }
    std::size_t find(std::size_t i) {
        if (parent_[i] != i) parent_[i] = find(parent_[i]);
        return parent_[i];
    }
    void join(std::size_t i, std::size_t j) {
        i = find(i);
        j = find(j);
        if (i != j) parent_[std::max(i, j)] = std::min(i, j);
    }
private:
    std::vector<std::size_t> parent_;
};

std::array<long double, 3> outward_area_normal(
    const TetrahedralMesh3D& mesh, const Key& key, std::size_t opposite
) {
    const auto& a = mesh.nodes[key[0]];
    const auto& b = mesh.nodes[key[1]];
    const auto& c = mesh.nodes[key[2]];
    const auto& o = mesh.nodes[opposite];
    const std::array<long double, 3> u{
        static_cast<long double>(b.x) - a.x,
        static_cast<long double>(b.y) - a.y,
        static_cast<long double>(b.z) - a.z
    };
    const std::array<long double, 3> v{
        static_cast<long double>(c.x) - a.x,
        static_cast<long double>(c.y) - a.y,
        static_cast<long double>(c.z) - a.z
    };
    std::array<long double, 3> normal{
        (u[1]*v[2] - u[2]*v[1])/2.0L,
        (u[2]*v[0] - u[0]*v[2])/2.0L,
        (u[0]*v[1] - u[1]*v[0])/2.0L
    };
    const long double orientation =
        normal[0]*(static_cast<long double>(o.x)-a.x) +
        normal[1]*(static_cast<long double>(o.y)-a.y) +
        normal[2]*(static_cast<long double>(o.z)-a.z);
    if (!std::isfinite(orientation) || orientation == 0.0L) {
        throw std::overflow_error("Cannot orient conservative flux face.");
    }
    if (orientation > 0) for (auto& x : normal) x = -x;
    return normal;
}

} // namespace

ConservativeFaceFlux3DResult reconstruct_conservative_face_flux_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const P1TetrahedralFlux3D> raw_cell_flux,
    std::span<const double> cell_source_heat_rate,
    std::span<const PrescribedBoundaryHeatRate3D> prescribed_boundary_rates,
    std::span<const BoundaryFace3D> adjustable_boundary_faces
) {
    mesh.validate();
    const std::size_t n = mesh.elements.size();
    if (raw_cell_flux.size() != n || cell_source_heat_rate.size() != n) {
        throw std::invalid_argument("Flux reconstruction needs one raw flux and source per tetrahedron.");
    }
    for (std::size_t i = 0; i < n; ++i) {
        const auto& q = raw_cell_flux[i].flux;
        if (!std::isfinite(q.x) || !std::isfinite(q.y) || !std::isfinite(q.z) ||
            !std::isfinite(cell_source_heat_rate[i])) {
            throw std::invalid_argument("Nonfinite raw flux or cell heat source.");
        }
    }

    std::map<Key, Face> faces;
    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto& tet = mesh.elements[cell];
        for (std::size_t omit = 0; omit < 4; ++omit) {
            Key face{};
            std::size_t k = 0;
            for (std::size_t j = 0; j < 4; ++j) {
                if (j != omit) face[k++] = tet[j];
            }
            face = sorted(face);
            auto [it, inserted] = faces.try_emplace(face, Face{cell, tet[omit]});
            if (!inserted) {
                if (it->second.neighbor != no_cell) {
                    throw std::invalid_argument("Nonmanifold tetrahedral face.");
                }
                it->second.neighbor = cell;
            }
        }
    }

    for (auto& [key, face] : faces) {
        const auto normal = outward_area_normal(mesh, key, face.opposite);
        face.area = checked(std::hypot(normal[0], normal[1], normal[2]),
                            "Face area overflow.");
        if (!(face.area > 0.0)) {
            throw std::invalid_argument("Conservative face area must be positive.");
        }
        if (face.neighbor != no_cell) {
            // Average one-sided physical heat-flux vectors into a single
            // orientation-independent interior interface heat rate.
            const auto& q0 = raw_cell_flux[face.cell].flux;
            const auto& q1 = raw_cell_flux[face.neighbor].flux;
            face.rate = checked(
                (normal[0]*(static_cast<long double>(q0.x)+q1.x) +
                 normal[1]*(static_cast<long double>(q0.y)+q1.y) +
                 normal[2]*(static_cast<long double>(q0.z)+q1.z))/2.0L,
                "Interior face heat rate overflow.");
        }
    }

    std::set<Key> claimed;
    for (const auto& bc : prescribed_boundary_rates) {
        const auto key = sorted(bc.face);
        const auto it = faces.find(key);
        if (it == faces.end() || it->second.neighbor != no_cell ||
            !claimed.insert(key).second) {
            throw std::invalid_argument("Prescribed heat rate must be on a unique exterior face.");
        }
        if (!std::isfinite(bc.outward_heat_rate)) {
            throw std::invalid_argument("Prescribed boundary heat rate must be finite.");
        }
        it->second.rate = bc.outward_heat_rate;
    }
    for (const auto& unassigned : adjustable_boundary_faces) {
        const auto key = sorted(unassigned);
        const auto it = faces.find(key);
        if (it == faces.end() || it->second.neighbor != no_cell ||
            !claimed.insert(key).second) {
            throw std::invalid_argument("Adjustable face must be a unique unused exterior face.");
        }
        const auto normal = outward_area_normal(mesh, key, it->second.opposite);
        const auto& q = raw_cell_flux[it->second.cell].flux;
        it->second.rate = checked(normal[0]*q.x + normal[1]*q.y + normal[2]*q.z,
                                  "Initial Dirichlet face heat rate overflow.");
        it->second.adjustable = true;
    }

    DisjointSet components(n);
    for (const auto& [key, face] : faces) {
        static_cast<void>(key);
        if (face.neighbor != no_cell) components.join(face.cell, face.neighbor);
    }
    std::vector<long double> rhs(n, 0.0L);
    std::vector<long double> component_defect(n, 0.0L);
    std::vector<std::size_t> representative(n, no_cell);
    std::vector<bool> has_open_boundary(n, false);
    std::vector<long double> component_scale(n, 1.0L);
    for (std::size_t cell = 0; cell < n; ++cell) {
        rhs[cell] = cell_source_heat_rate[cell];
        component_scale[components.find(cell)] += std::abs(rhs[cell]);
    }
    for (const auto& [key, face] : faces) {
        static_cast<void>(key);
        rhs[face.cell] -= face.rate;
        if (face.neighbor != no_cell) rhs[face.neighbor] += face.rate;
        if (face.adjustable) has_open_boundary[components.find(face.cell)] = true;
        if (face.neighbor == no_cell && !face.adjustable) {
            component_scale[components.find(face.cell)] +=
                std::abs(static_cast<long double>(face.rate));
        }
    }
    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto root = components.find(cell);
        representative[root] = std::min(representative[root], cell);
        component_defect[root] += rhs[cell];
    }
    std::vector<bool> pinned(n, false);
    for (std::size_t root = 0; root < n; ++root) {
        if (representative[root] == no_cell || has_open_boundary[root]) continue;
        if (std::abs(component_defect[root]) > 1e-11L * component_scale[root]) {
            throw std::invalid_argument(
                "Unanchored component has an incompatible total source/flux budget.");
        }
        pinned[representative[root]] = true; // numerical graph gauge, not heat flow
    }

    // Minimize sum_f (delta Q_f)^2 / area_f subject to B*Q=source.
    // L = B diag(area) B^T, with one gauge pin on closed components.
    std::vector<linalg::Triplet> entries;
    entries.reserve(faces.size()*4 + n);
    for (const auto& [key, face] : faces) {
        static_cast<void>(key);
        if (face.neighbor == no_cell && !face.adjustable) continue;
        if (!pinned[face.cell]) {
            entries.push_back({face.cell, face.cell, face.area});
        }
        if (face.neighbor != no_cell) {
            if (!pinned[face.neighbor]) {
                entries.push_back({face.neighbor, face.neighbor, face.area});
            }
            if (!pinned[face.cell] && !pinned[face.neighbor]) {
                entries.push_back({face.cell, face.neighbor, -face.area});
                entries.push_back({face.neighbor, face.cell, -face.area});
            }
        }
    }
    std::vector<double> right(n, 0.0);
    for (std::size_t cell = 0; cell < n; ++cell) {
        if (pinned[cell]) entries.push_back({cell, cell, 1.0});
        else right[cell] = checked(rhs[cell], "Face conservation rhs overflow.");
    }
    const auto matrix = linalg::CsrMatrix::from_triplets(n, n, entries);
    linalg::SolverOptions options;
    options.max_iterations = std::max<std::size_t>(
        2000, std::min<std::size_t>(n, 20000)*8);
    options.relative_tolerance = 1e-12;
    options.absolute_tolerance = 1e-14;
    options.preconditioner = linalg::PreconditionerType::Jacobi;
    const auto solver = linalg::solve_conjugate_gradient(matrix, right, options);
    if (!solver.converged) {
        throw std::runtime_error("Sparse face conservation correction failed: " + solver.message);
    }
    for (auto& [key, face] : faces) {
        static_cast<void>(key);
        if (face.neighbor == no_cell && !face.adjustable) continue;
        long double difference = solver.solution[face.cell];
        if (face.neighbor != no_cell) difference -= solver.solution[face.neighbor];
        face.rate = checked(
            static_cast<long double>(face.rate) + face.area*difference,
            "Corrected face heat rate overflow.");
    }

    ConservativeFaceFlux3DResult out;
    out.faces.reserve(faces.size());
    out.cell_balance_defect.resize(n);
    std::vector<long double> defect(n, 0.0L);
    std::vector<long double> local_scale(n, 1.0L);
    for (std::size_t cell = 0; cell < n; ++cell) {
        defect[cell] = cell_source_heat_rate[cell];
        local_scale[cell] += std::abs(static_cast<long double>(cell_source_heat_rate[cell]));
    }
    for (const auto& [key, face] : faces) {
        defect[face.cell] -= face.rate;
        if (face.neighbor != no_cell) defect[face.neighbor] += face.rate;
        local_scale[face.cell] += std::abs(static_cast<long double>(face.rate));
        if (face.neighbor != no_cell) {
            local_scale[face.neighbor] += std::abs(static_cast<long double>(face.rate));
        }
        out.faces.push_back({key, face.cell,
            face.neighbor == no_cell ? std::nullopt : std::optional<std::size_t>{face.neighbor},
            face.area, face.rate});
    }
    for (std::size_t cell = 0; cell < n; ++cell) {
        out.cell_balance_defect[cell] = checked(defect[cell], "Cell heat defect overflow.");
        out.maximum_absolute_cell_defect = std::max(
            out.maximum_absolute_cell_defect, std::abs(out.cell_balance_defect[cell]));
    }
    for (std::size_t cell = 0; cell < n; ++cell) {
        if (std::abs(defect[cell]) > 1e-9L * local_scale[cell]) {
            throw std::runtime_error(
                "Face equilibration did not converge to per-cell conservation.");
        }
    }
    out.correction_iterations = solver.iterations;
    return out;
}

ConservativeFaceFlux3DResult reconstruct_conservative_face_flux_3d(
    const ScalarDiffusion3DProblem& problem, std::span<const double> nodal_values
) {
    const auto& mesh = problem.mesh;
    const auto raw = evaluate_p1_tetrahedral_flux_3d(problem, nodal_values);
    if (!problem.source ||
        (!problem.cell_source.empty() && problem.cell_source.size() != mesh.elements.size())) {
        throw std::invalid_argument("Reconstruction requires valid spatial and per-cell sources.");
    }
    std::vector<double> cell_source(mesh.elements.size(), 0.0);
    for (std::size_t cell = 0; cell < mesh.elements.size(); ++cell) {
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, cell);
        const auto& tet = mesh.elements[cell];
        const double jacobian = std::abs(basis.signed_jacobian_determinant);
        const double extra = problem.cell_source.empty() ? 0.0 : problem.cell_source[cell];
        if (!std::isfinite(extra)) throw std::invalid_argument("Nonfinite cell source.");
        long double rate = 0.0L;
        for (const auto& q : tetrahedron_degree_two_quadrature) {
            const auto pt = detail::simplex_point(mesh, tet, q.shape);
            const double value = problem.source(pt.x, pt.y, pt.z) + extra;
            if (!std::isfinite(value)) {
                throw std::invalid_argument("Nonfinite volumetric heat source.");
            }
            rate += static_cast<long double>(value)*jacobian*q.weight;
        }
        cell_source[cell] = checked(rate, "Integrated cell generation overflow.");
    }

    detail::ExteriorFaceAssignments3D claims(mesh);
    std::vector<PrescribedBoundaryHeatRate3D> prescribed;
    std::vector<BoundaryFace3D> adjustable;
    std::vector<bool> fixed_nodes(mesh.nodes.size(), false);
    std::set<Key> assigned;
    for (const auto& fixed : problem.dirichlet_constraints) {
        if (fixed.node >= mesh.nodes.size() || !std::isfinite(fixed.value)) {
            throw std::invalid_argument("Invalid Dirichlet node constraint.");
        }
        fixed_nodes[fixed.node] = true;
    }
    for (const auto& bc : problem.dirichlet_faces) {
        claims.claim(bc.face);
        if (!bc.value) throw std::invalid_argument("Missing Dirichlet face function.");
        assigned.insert(sorted(bc.face));
        adjustable.push_back(bc.face);
        for (auto node : bc.face) fixed_nodes[node] = true;
    }
    const auto integrate = [&](const BoundaryFace3D& face, const auto& sample) {
        const double jacobian = detail::face_jacobian(mesh, face);
        long double total = 0.0L;
        for (const auto& q : triangle_face_degree_two_quadrature) {
            const auto pt = detail::simplex_point(mesh, face, q.shape);
            total += static_cast<long double>(sample(pt, q.shape))*jacobian*q.weight;
        }
        return checked(total, "Boundary heat rate overflow.");
    };
    for (const auto& bc : problem.neumann_faces) {
        claims.claim(bc.face);
        if (!bc.outward_flux) throw std::invalid_argument("Missing Neumann flux.");
        assigned.insert(sorted(bc.face));
        prescribed.push_back({bc.face, integrate(bc.face, [&](const Point3D& pt, const auto&) {
            const double q = bc.outward_flux(pt.x, pt.y, pt.z);
            if (!std::isfinite(q)) throw std::invalid_argument("Nonfinite Neumann flux.");
            return q;
        })});
    }
    for (const auto& bc : problem.robin_faces) {
        claims.claim(bc.face);
        if (!bc.transfer_coefficient || !bc.ambient_value) {
            throw std::invalid_argument("Missing Robin boundary function.");
        }
        assigned.insert(sorted(bc.face));
        prescribed.push_back({bc.face, integrate(bc.face, [&](const Point3D& pt, const auto& shape) {
            const double h = bc.transfer_coefficient(pt.x, pt.y, pt.z);
            const double ambient = bc.ambient_value(pt.x, pt.y, pt.z);
            if (!std::isfinite(h) || h < 0.0 || !std::isfinite(ambient)) {
                throw std::invalid_argument("Invalid Robin transfer or ambient temperature.");
            }
            long double temp = 0.0L;
            for (std::size_t j = 0; j < 3; ++j) {
                temp += static_cast<long double>(nodal_values[bc.face[j]])*shape[j];
            }
            return checked(h*(temp-ambient), "Robin heat-flux density overflow.");
        })});
    }
    for (const auto& face : extract_boundary_faces_3d(mesh)) {
        if (assigned.contains(sorted(face))) continue;
        if (fixed_nodes[face[0]] && fixed_nodes[face[1]] && fixed_nodes[face[2]]) {
            adjustable.push_back(face);
        }
    }
    return reconstruct_conservative_face_flux_3d(mesh, raw, cell_source,
                                                  prescribed, adjustable);
}

} // namespace nova::fem
