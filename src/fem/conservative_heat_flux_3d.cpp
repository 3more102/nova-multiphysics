#include "nova/fem/conservative_heat_flux_3d.hpp"

#include "nova/fem/energy_balance_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"
#include "nova/fem/simplex_quadrature.hpp"
#include "nova/linalg/conjugate_gradient.hpp"
#include "p1_integration_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace nova::fem {
namespace {

enum class FaceKind { Interior, Insulated, Dirichlet, Neumann, Robin };

struct FaceWork {
    ConservativeFaceHeat3D data;
    std::size_t opposite{};
    FaceKind kind{FaceKind::Insulated};
    double raw_rate{};
};

struct DisjointSet {
    explicit DisjointSet(std::size_t n) : parents(n) {
        std::iota(parents.begin(), parents.end(), std::size_t{0});
    }

    std::size_t find(std::size_t x) {
        while (parents[x] != x) {
            parents[x] = parents[parents[x]];
            x = parents[x];
        }
        return x;
    }

    void join(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a != b) parents[b] = a;
    }

    std::vector<std::size_t> parents;
};

[[nodiscard]] BoundaryFace3D key(BoundaryFace3D face) {
    std::sort(face.begin(), face.end());
    return face;
}

[[nodiscard]] double checked(long double value, const char* context) {
    if (!std::isfinite(value) ||
        std::abs(value) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error(context);
    }
    return static_cast<double>(value);
}

[[nodiscard]] std::array<long double, 3> outward_area_vector(
    const TetrahedralMesh3D& mesh, const BoundaryFace3D& face,
    std::size_t opposite
) {
    const auto& a = mesh.nodes[face[0]];
    const auto& b = mesh.nodes[face[1]];
    const auto& c = mesh.nodes[face[2]];
    const auto& d = mesh.nodes[opposite];
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
    std::array<long double, 3> area{
        (u[1]*v[2] - u[2]*v[1]) / 2,
        (u[2]*v[0] - u[0]*v[2]) / 2,
        (u[0]*v[1] - u[1]*v[0]) / 2
    };
    const long double into_cell =
        area[0]*(static_cast<long double>(d.x)-a.x) +
        area[1]*(static_cast<long double>(d.y)-a.y) +
        area[2]*(static_cast<long double>(d.z)-a.z);
    if (!std::isfinite(into_cell) || into_cell == 0.0L) {
        throw std::overflow_error("Conservative heat face normal is not representable.");
    }
    if (into_cell > 0) {
        for (auto& v_component : area) v_component = -v_component;
    }
    return area;
}

[[nodiscard]] double flux_rate(
    const std::array<long double, 3>& area, const Point3D& flux
) {
    return checked(area[0]*flux.x + area[1]*flux.y + area[2]*flux.z,
                   "Raw conductive face heat rate overflowed.");
}

} // namespace

ConservativeHeatFluxResult3D reconstruct_conservative_heat_flux_3d(
    const ScalarDiffusion3DProblem& problem,
    const std::span<const double> nodal_values
) {
    // Reuse the established weak-form audit's complete validation of fields,
    // material selection, sources, exterior face exclusivity and fixed nodes.
    (void)evaluate_diffusion_energy_balance_3d(problem, nodal_values);
    const auto& mesh = problem.mesh;
    const std::size_t n = mesh.elements.size();
    const auto cell_flux = evaluate_p1_tetrahedral_flux_3d(problem, nodal_values);

    // Create each unique face once. mesh.validate() (above) guarantees at most
    // two opposite-sided tetrahedra per face and valid vertex references.
    std::map<BoundaryFace3D, FaceWork> unique_faces;
    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto& tetra = mesh.elements[cell];
        for (std::size_t omitted = 0; omitted < 4; ++omitted) {
            BoundaryFace3D face{};
            std::size_t at = 0;
            for (std::size_t j = 0; j < 4; ++j) {
                if (j != omitted) face[at++] = tetra[j];
            }
            const auto canonical = key(face);
            auto [it, inserted] = unique_faces.try_emplace(canonical);
            if (inserted) {
                it->second.data.vertices = canonical;
                it->second.data.owner_cell = cell;
                it->second.opposite = tetra[omitted];
            } else {
                it->second.data.neighbor_cell = cell;
                it->second.kind = FaceKind::Interior;
            }
        }
    }

    std::map<BoundaryFace3D, FaceKind> boundaries;
    for (const auto& bc : problem.dirichlet_faces) {
        boundaries.emplace(key(bc.face), FaceKind::Dirichlet);
    }
    for (const auto& bc : problem.neumann_faces) {
        boundaries.emplace(key(bc.face), FaceKind::Neumann);
    }
    for (const auto& bc : problem.robin_faces) {
        boundaries.emplace(key(bc.face), FaceKind::Robin);
    }

    std::map<BoundaryFace3D, double> imposed;
    const auto integrate_face = [&](const BoundaryFace3D& face, auto&& density) {
        const double jacobian = detail::face_jacobian(mesh, face);
        long double total = 0.0L;
        for (const auto& q : triangle_face_degree_two_quadrature) {
            const auto p = detail::simplex_point(mesh, face, q.shape);
            total += static_cast<long double>(density(p, q.shape)) *
                     (jacobian * q.weight);
        }
        return checked(total, "Conservative heat boundary integral overflowed.");
    };
    for (const auto& bc : problem.neumann_faces) {
        imposed.emplace(key(bc.face), integrate_face(bc.face,
            [&](const Point3D& p, const auto&) {
                const double value = bc.outward_flux(p.x, p.y, p.z);
                if (!std::isfinite(value)) {
                    throw std::invalid_argument("Nonfinite Neumann flux.");
                }
                return value;
            }));
    }
    for (const auto& bc : problem.robin_faces) {
        imposed.emplace(key(bc.face), integrate_face(bc.face,
            [&](const Point3D& p, const auto& shape) {
                const double h = bc.transfer_coefficient(p.x, p.y, p.z);
                const double ambient = bc.ambient_value(p.x, p.y, p.z);
                if (!std::isfinite(h) || h < 0 || !std::isfinite(ambient)) {
                    throw std::invalid_argument("Invalid Robin boundary data.");
                }
                long double temp = 0.0L;
                for (std::size_t j = 0; j < 3; ++j) {
                    temp += static_cast<long double>(nodal_values[bc.face[j]]) * shape[j];
                }
                return checked(static_cast<long double>(h)*(temp - ambient),
                               "Robin heat density overflowed.");
            }));
    }

    ConservativeHeatFluxResult3D result;
    result.generated_heat_per_cell.assign(n, 0.0);
    result.outward_heat_per_cell.assign(n, 0.0);
    std::vector<long double> rhs(n, 0.0L);
    std::vector<long double> magnitude(n, 0.0L);

    // Use the same four-point volume quadrature and additive cell source as
    // the steady FEM assembly, including sign-changing and anisotropic cases.
    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto& tetra = mesh.elements[cell];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, cell);
        const double jacobian = std::abs(basis.signed_jacobian_determinant);
        const double regional = problem.cell_source.empty() ? 0 : problem.cell_source[cell];
        long double total = 0.0L;
        for (const auto& q : tetrahedron_degree_two_quadrature) {
            const auto p = detail::simplex_point(mesh, tetra, q.shape);
            const double source = problem.source(p.x, p.y, p.z) + regional;
            if (!std::isfinite(source)) {
                throw std::invalid_argument("Nonfinite total heat source.");
            }
            total += static_cast<long double>(source) * (jacobian * q.weight);
        }
        result.generated_heat_per_cell[cell] =
            checked(total, "Cell heat generation overflowed.");
        rhs[cell] = total;
        magnitude[cell] += std::abs(total);
    }

    DisjointSet components(n);
    for (auto& [face, work] : unique_faces) {
        const auto& info = work.data;
        const auto area = outward_area_vector(mesh, face, work.opposite);
        const double owner_rate = flux_rate(area, cell_flux[info.owner_cell].flux);
        if (info.neighbor_cell != ConservativeFaceHeat3D::no_neighbor) {
            // Owner and neighbor outward normals are opposites.
            const double neighbor_rate = flux_rate(area, cell_flux[info.neighbor_cell].flux);
            work.raw_rate = checked(
                (static_cast<long double>(owner_rate) + neighbor_rate) / 2,
                "Interior face target heat rate overflowed.");
            components.join(info.owner_cell, info.neighbor_cell);
        } else {
            const auto found = boundaries.find(face);
            work.kind = found == boundaries.end() ? FaceKind::Insulated : found->second;
            if (work.kind == FaceKind::Neumann || work.kind == FaceKind::Robin) {
                work.raw_rate = imposed.at(face);
            } else if (work.kind == FaceKind::Dirichlet) {
                work.raw_rate = owner_rate;
            } else {
                work.raw_rate = 0.0;
            }
        }
        rhs[info.owner_cell] -= work.raw_rate;
        magnitude[info.owner_cell] += std::abs(work.raw_rate);
        if (info.neighbor_cell != ConservativeFaceHeat3D::no_neighbor) {
            rhs[info.neighbor_cell] += work.raw_rate;
            magnitude[info.neighbor_cell] += std::abs(work.raw_rate);
        }
    }

    // The only adjustable boundary outlets are explicitly tagged Dirichlet
    // faces. A component without one needs exact prescribed net balance.
    std::vector<bool> anchored(n, false), pinned(n, false);
    for (const auto& [face, work] : unique_faces) {
        (void)face;
        if (work.kind == FaceKind::Dirichlet) {
            anchored[components.find(work.data.owner_cell)] = true;
        }
    }
    std::vector<long double> component_rhs(n, 0.0L), component_scale(n, 0.0L);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t root = components.find(i);
        component_rhs[root] += rhs[i];
        component_scale[root] += magnitude[i];
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (components.find(i) != i || anchored[i]) continue;
        const long double tolerance = 1e-10L *
            std::max(component_scale[i], 1e-12L);
        if (std::abs(component_rhs[i]) > tolerance) {
            throw std::invalid_argument(
                "Heat source and prescribed face flux are incompatible without a Dirichlet face outlet.");
        }
        // Fix exactly one potential per unanchored component (gauge only).
        pinned[i] = true;
    }

    std::vector<double> diagonal(n, 0.0);
    std::vector<linalg::Triplet> entries;
    entries.reserve(n * 5);
    for (const auto& [face, work] : unique_faces) {
        (void)face;
        const auto a = work.data.owner_cell;
        if (work.kind == FaceKind::Interior) {
            const auto b = work.data.neighbor_cell;
            diagonal[a] += 1;
            diagonal[b] += 1;
            if (!pinned[a] && !pinned[b]) {
                entries.push_back({a, b, -1});
                entries.push_back({b, a, -1});
            }
        } else if (work.kind == FaceKind::Dirichlet) {
            diagonal[a] += 1;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        entries.push_back({i, i, pinned[i] ? 1.0 : diagonal[i]});
        if (pinned[i]) rhs[i] = 0.0;
    }
    std::vector<double> right_hand_side;
    right_hand_side.reserve(n);
    for (long double value : rhs) {
        right_hand_side.push_back(checked(value, "Face-flow correction RHS overflowed."));
    }

    const auto matrix = linalg::CsrMatrix::from_triplets(n, n, entries);
    linalg::SolverOptions options;
    options.max_iterations = std::max<std::size_t>(
        1000, std::min<std::size_t>(20000, n > 3333 ? 20000 : 6*n)
    );
    options.relative_tolerance = 1e-12;
    options.absolute_tolerance = 1e-14;
    options.preconditioner = linalg::PreconditionerType::Jacobi;
    const auto correction = linalg::solve_conjugate_gradient(
        matrix, right_hand_side, options
    );
    if (!correction.converged) {
        throw std::runtime_error(
            "Conservative heat-flux correction did not converge: " + correction.message
        );
    }
    result.correction_iterations = correction.iterations;

    std::vector<long double> outward(n, 0.0L);
    result.faces.reserve(unique_faces.size());
    for (const auto& [face, work] : unique_faces) {
        (void)face;
        auto output = work.data;
        long double rate = work.raw_rate;
        if (work.kind == FaceKind::Interior) {
            rate += static_cast<long double>(correction.solution[output.owner_cell]) -
                    correction.solution[output.neighbor_cell];
        } else if (work.kind == FaceKind::Dirichlet) {
            rate += correction.solution[output.owner_cell];
        }
        output.outward_heat_rate = checked(rate, "Corrected face heat rate overflowed.");
        outward[output.owner_cell] += rate;
        if (output.neighbor_cell != ConservativeFaceHeat3D::no_neighbor) {
            outward[output.neighbor_cell] -= rate;
        }
        result.faces.push_back(output);
    }

    long double maximum = 0.0L, scale = 1.0L;
    for (std::size_t i = 0; i < n; ++i) {
        result.outward_heat_per_cell[i] = checked(outward[i], "Cell heat rate overflowed.");
        maximum = std::max(maximum, std::abs(
            static_cast<long double>(result.generated_heat_per_cell[i]) - outward[i]));
        scale = std::max(scale, magnitude[i]);
    }
    result.max_cell_balance_defect = checked(maximum, "Cell balance defect overflowed.");
    if (maximum > 1e-9L * scale) {
        throw std::runtime_error("Conservative face-flow reconstruction has excessive local residual.");
    }
    return result;
}

} // namespace nova::fem
