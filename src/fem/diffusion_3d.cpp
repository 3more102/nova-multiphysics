#include "nova/fem/diffusion_3d.hpp"

#include "nova/fem/tetrahedron_h1.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using FaceKey = BoundaryFace3D;

struct TetrahedronQuadraturePoint {
    double r{};
    double s{};
    double t{};
    double weight{};
};

struct FaceQuadraturePoint {
    std::array<double, 3> shape{};
    double weight{};
};

[[nodiscard]] FaceKey canonical_face(FaceKey face) {
    std::sort(face.begin(), face.end());
    return face;
}

// Positive, symmetric, degree-two exact rule on the reference tetrahedron.
// The four weights sum to its volume, 1/6.
[[nodiscard]] constexpr std::array<TetrahedronQuadraturePoint, 4>
tetrahedron_degree_two_quadrature() {
    constexpr double a = 0.5854101966249685;
    constexpr double b = 0.1381966011250105;
    constexpr double weight = 1.0 / 24.0;
    return {{
        {b, b, b, weight},
        {a, b, b, weight},
        {b, a, b, weight},
        {b, b, a, weight}
    }};
}

// Positive, symmetric, degree-two exact rule on a reference triangle.
// The weights sum to its area, 1/2.
[[nodiscard]] constexpr std::array<FaceQuadraturePoint, 3>
triangle_degree_two_quadrature() {
    constexpr double a = 2.0 / 3.0;
    constexpr double b = 1.0 / 6.0;
    constexpr double weight = 1.0 / 6.0;
    return {{
        {{a, b, b}, weight},
        {{b, a, b}, weight},
        {{b, b, a}, weight}
    }};
}

void prescribe(
    std::vector<std::optional<double>>& prescribed,
    const std::size_t node,
    const double value
) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("Dirichlet values must be finite.");
    }
    if (prescribed[node].has_value()) {
        const double previous = *prescribed[node];
        const double scale = std::max({
            1.0, std::abs(previous), std::abs(value)
        });
        if (std::abs(previous - value) >
            64.0 * std::numeric_limits<double>::epsilon() * scale) {
            throw std::invalid_argument(
                "Conflicting Dirichlet values at a shared 3D vertex."
            );
        }
        return;
    }
    prescribed[node] = value;
}

[[nodiscard]] double twice_face_area(
    const Point3D& a,
    const Point3D& b,
    const Point3D& c
) {
    const double ux = b.x - a.x;
    const double uy = b.y - a.y;
    const double uz = b.z - a.z;
    const double vx = c.x - a.x;
    const double vy = c.y - a.y;
    const double vz = c.z - a.z;
    const double nx = uy * vz - uz * vy;
    const double ny = uz * vx - ux * vz;
    const double nz = ux * vy - uy * vx;
    return std::hypot(std::hypot(nx, ny), nz);
}

[[nodiscard]] Point3D interpolate_face_point(
    const TetrahedralMesh3D& mesh,
    const BoundaryFace3D& face,
    const std::array<double, 3>& shape
) {
    const Point3D& a = mesh.nodes[face[0]];
    const Point3D& b = mesh.nodes[face[1]];
    const Point3D& c = mesh.nodes[face[2]];
    return {
        shape[0] * a.x + shape[1] * b.x + shape[2] * c.x,
        shape[0] * a.y + shape[1] * b.y + shape[2] * c.y,
        shape[0] * a.z + shape[1] * b.z + shape[2] * c.z
    };
}

} // namespace

ScalarDiffusion3DResult solve_scalar_diffusion_3d(
    const ScalarDiffusion3DProblem& problem,
    const linalg::SolverOptions& solver_options
) {
    problem.mesh.validate();
    const H1FunctionSpace3D space(problem.mesh);
    const std::size_t node_count = space.dimension();

    if (!problem.diffusivity || !problem.source) {
        throw std::invalid_argument(
            "3D diffusivity and source functions are required."
        );
    }
    if (problem.dirichlet_constraints.empty() &&
        problem.dirichlet_faces.empty()) {
        throw std::invalid_argument(
            "At least one 3D Dirichlet constraint or boundary face is required."
        );
    }

    std::set<FaceKey> exterior_faces;
    for (const BoundaryFace3D& face : extract_boundary_faces_3d(problem.mesh)) {
        exterior_faces.insert(canonical_face(face));
    }

    const auto validate_exterior_face = [&](const BoundaryFace3D& face) {
        for (const std::size_t node : face) {
            if (node >= node_count) {
                throw std::out_of_range(
                    "3D boundary face references a missing mesh vertex."
                );
            }
        }
        if (face[0] == face[1] || face[1] == face[2] ||
            face[2] == face[0]) {
            throw std::invalid_argument(
                "A 3D boundary face must contain three distinct vertices."
            );
        }
        const FaceKey key = canonical_face(face);
        if (exterior_faces.count(key) == 0U) {
            throw std::invalid_argument(
                "3D boundary conditions require an exterior tetrahedral face."
            );
        }
        return key;
    };

    std::vector<std::optional<double>> prescribed(node_count);
    std::set<std::size_t> explicit_nodes;
    for (const DirichletConstraint3D& constraint :
         problem.dirichlet_constraints) {
        if (constraint.node >= node_count) {
            throw std::out_of_range(
                "3D Dirichlet constraint references a missing vertex."
            );
        }
        if (!explicit_nodes.insert(constraint.node).second) {
            throw std::invalid_argument(
                "Duplicate explicit 3D Dirichlet vertex constraint."
            );
        }
        prescribe(prescribed, constraint.node, constraint.value);
    }

    std::set<FaceKey> dirichlet_faces;
    for (const DirichletBoundaryFace3D& condition :
         problem.dirichlet_faces) {
        const FaceKey key = validate_exterior_face(condition.face);
        if (!condition.value) {
            throw std::invalid_argument(
                "3D Dirichlet face requires a value function."
            );
        }
        if (!dirichlet_faces.insert(key).second) {
            throw std::invalid_argument("Duplicate 3D Dirichlet boundary face.");
        }
        for (const std::size_t node : condition.face) {
            const Point3D& point = problem.mesh.nodes[node];
            prescribe(
                prescribed, node,
                condition.value(point.x, point.y, point.z)
            );
        }
    }

    std::set<FaceKey> natural_faces;
    for (const NeumannBoundaryFace3D& condition : problem.neumann_faces) {
        const FaceKey key = validate_exterior_face(condition.face);
        if (!condition.outward_flux) {
            throw std::invalid_argument(
                "3D Neumann face requires an outward-flux function."
            );
        }
        if (dirichlet_faces.count(key) != 0U ||
            !natural_faces.insert(key).second) {
            throw std::invalid_argument(
                "Duplicate or incompatible 3D boundary-face conditions."
            );
        }
    }

    for (const RobinBoundaryFace3D& condition : problem.robin_faces) {
        const FaceKey key = validate_exterior_face(condition.face);
        if (!condition.transfer_coefficient || !condition.ambient_value) {
            throw std::invalid_argument(
                "3D Robin face requires transfer and ambient functions."
            );
        }
        if (dirichlet_faces.count(key) != 0U ||
            !natural_faces.insert(key).second) {
            throw std::invalid_argument(
                "Duplicate or incompatible 3D boundary-face conditions."
            );
        }
    }

    constexpr std::size_t local_dof_count = 4U;
    constexpr std::size_t volume_entries_per_element =
        local_dof_count * local_dof_count;
    constexpr std::size_t robin_entries_per_face =
        3U * local_dof_count - 3U; // overwritten below with exact face count
    static_cast<void>(robin_entries_per_face);
    constexpr std::size_t robin_matrix_entries_per_face = 3U * 3U * 3U;

    std::vector<linalg::Triplet> entries;
    const std::size_t capacity = entries.max_size();
    if (problem.mesh.elements.size() >
        capacity / volume_entries_per_element) {
        throw std::length_error("3D finite-element volume matrix is too large.");
    }
    const std::size_t volume_entry_count =
        problem.mesh.elements.size() * volume_entries_per_element;
    if (problem.robin_faces.size() >
        (capacity - volume_entry_count) / robin_matrix_entries_per_face) {
        throw std::length_error("3D Robin boundary matrix is too large.");
    }
    entries.reserve(
        volume_entry_count +
        problem.robin_faces.size() * robin_matrix_entries_per_face
    );

    std::vector<double> rhs(node_count, 0.0);
    constexpr auto volume_quadrature = tetrahedron_degree_two_quadrature();

    for (std::size_t element = 0;
         element < problem.mesh.elements.size();
         ++element) {
        const Tetrahedron3D& cell = problem.mesh.elements[element];
        const Point3D& p0 = problem.mesh.nodes[cell[0]];
        const Point3D& p1 = problem.mesh.nodes[cell[1]];
        const Point3D& p2 = problem.mesh.nodes[cell[2]];
        const Point3D& p3 = problem.mesh.nodes[cell[3]];

        const double j00 = p1.x - p0.x;
        const double j01 = p2.x - p0.x;
        const double j02 = p3.x - p0.x;
        const double j10 = p1.y - p0.y;
        const double j11 = p2.y - p0.y;
        const double j12 = p3.y - p0.y;
        const double j20 = p1.z - p0.z;
        const double j21 = p2.z - p0.z;
        const double j22 = p3.z - p0.z;
        const auto physical = evaluate_tetrahedron_p1_physical_gradients_3d(
            problem.mesh, element
        );
        const double absolute_determinant =
            std::abs(physical.signed_jacobian_determinant);

        std::array<std::array<double, 4>, 4> local_matrix{};
        std::array<double, 4> local_rhs{};

        for (const TetrahedronQuadraturePoint& quadrature :
             volume_quadrature) {
            const auto basis = evaluate_tetrahedron_p1_basis_3d(
                quadrature.r, quadrature.s, quadrature.t
            );
            const double x = p0.x +
                j00 * quadrature.r + j01 * quadrature.s +
                j02 * quadrature.t;
            const double y = p0.y +
                j10 * quadrature.r + j11 * quadrature.s +
                j12 * quadrature.t;
            const double z = p0.z +
                j20 * quadrature.r + j21 * quadrature.s +
                j22 * quadrature.t;

            const double diffusivity = problem.diffusivity(x, y, z);
            const double source = problem.source(x, y, z);
            if (!std::isfinite(diffusivity) || diffusivity <= 0.0) {
                throw std::invalid_argument(
                    "3D diffusivity must be finite and strictly positive."
                );
            }
            if (!std::isfinite(source)) {
                throw std::invalid_argument("3D source must be finite.");
            }

            const double weight = absolute_determinant * quadrature.weight;
            for (std::size_t i = 0; i < local_dof_count; ++i) {
                local_rhs[i] += source * basis.values[i] * weight;
                for (std::size_t j = 0; j < local_dof_count; ++j) {
                    const double gradient_dot =
                        physical.gradients[i][0] * physical.gradients[j][0] +
                        physical.gradients[i][1] * physical.gradients[j][1] +
                        physical.gradients[i][2] * physical.gradients[j][2];
                    local_matrix[i][j] +=
                        diffusivity * gradient_dot * weight;
                }
            }
        }

        for (std::size_t i = 0; i < local_dof_count; ++i) {
            const std::size_t global_i = cell[i];
            rhs[global_i] += local_rhs[i];
            if (!std::isfinite(rhs[global_i])) {
                throw std::overflow_error(
                    "3D element load-vector accumulation overflowed."
                );
            }
            for (std::size_t j = 0; j < local_dof_count; ++j) {
                if (!std::isfinite(local_matrix[i][j])) {
                    throw std::overflow_error(
                        "3D element stiffness entry is non-finite."
                    );
                }
                entries.push_back({
                    global_i, cell[j], local_matrix[i][j]
                });
            }
        }
    }

    constexpr auto face_quadrature = triangle_degree_two_quadrature();

    const auto integrate_face = [&](const BoundaryFace3D& face,
                                    const auto& callback,
                                    const bool robin) {
        const Point3D& a = problem.mesh.nodes[face[0]];
        const Point3D& b = problem.mesh.nodes[face[1]];
        const Point3D& c = problem.mesh.nodes[face[2]];
        const double area_scale = twice_face_area(a, b, c);
        if (!std::isfinite(area_scale) || area_scale <= 0.0) {
            throw std::invalid_argument(
                "3D boundary face area must be finite and positive."
            );
        }

        for (const FaceQuadraturePoint& quadrature : face_quadrature) {
            const Point3D point = interpolate_face_point(
                problem.mesh, face, quadrature.shape
            );
            const double weight = area_scale * quadrature.weight;
            if (robin) {
                const auto [transfer, ambient] = callback(
                    point.x, point.y, point.z
                );
                if (!std::isfinite(transfer) || transfer < 0.0) {
                    throw std::invalid_argument(
                        "3D Robin transfer coefficient must be finite and nonnegative."
                    );
                }
                if (!std::isfinite(ambient)) {
                    throw std::invalid_argument(
                        "3D Robin ambient value must be finite."
                    );
                }
                const double ambient_load = transfer * ambient;
                if (!std::isfinite(ambient_load)) {
                    throw std::overflow_error(
                        "3D Robin ambient load is non-finite."
                    );
                }
                for (std::size_t i = 0; i < 3U; ++i) {
                    const std::size_t global_i = face[i];
                    rhs[global_i] +=
                        ambient_load * quadrature.shape[i] * weight;
                    if (!std::isfinite(rhs[global_i])) {
                        throw std::overflow_error(
                            "3D Robin load-vector accumulation overflowed."
                        );
                    }
                    for (std::size_t j = 0; j < 3U; ++j) {
                        const double value =
                            transfer * quadrature.shape[i] *
                            quadrature.shape[j] * weight;
                        if (!std::isfinite(value)) {
                            throw std::overflow_error(
                                "3D Robin matrix entry is non-finite."
                            );
                        }
                        entries.push_back({global_i, face[j], value});
                    }
                }
            } else {
                const double flux = callback(point.x, point.y, point.z);
                if (!std::isfinite(flux)) {
                    throw std::invalid_argument(
                        "3D Neumann outward flux must be finite."
                    );
                }
                for (std::size_t i = 0; i < 3U; ++i) {
                    const std::size_t global_i = face[i];
                    rhs[global_i] -=
                        flux * quadrature.shape[i] * weight;
                    if (!std::isfinite(rhs[global_i])) {
                        throw std::overflow_error(
                            "3D Neumann load-vector accumulation overflowed."
                        );
                    }
                }
            }
        }
    };

    for (const NeumannBoundaryFace3D& condition : problem.neumann_faces) {
        integrate_face(
            condition.face,
            [&](const double x, const double y, const double z) {
                return condition.outward_flux(x, y, z);
            },
            false
        );
    }

    for (const RobinBoundaryFace3D& condition : problem.robin_faces) {
        integrate_face(
            condition.face,
            [&](const double x, const double y, const double z) {
                return std::pair{
                    condition.transfer_coefficient(x, y, z),
                    condition.ambient_value(x, y, z)
                };
            },
            true
        );
    }

    // Symmetric strong Dirichlet elimination preserves CG's SPD requirement.
    std::vector<linalg::Triplet> constrained_entries;
    if (node_count > constrained_entries.max_size() - entries.size()) {
        throw std::length_error("Constrained 3D finite-element matrix is too large.");
    }
    constrained_entries.reserve(entries.size() + node_count);
    for (const linalg::Triplet& entry : entries) {
        if (!prescribed[entry.row].has_value()) {
            if (prescribed[entry.column].has_value()) {
                rhs[entry.row] -=
                    entry.value * *prescribed[entry.column];
                if (!std::isfinite(rhs[entry.row])) {
                    throw std::overflow_error(
                        "3D Dirichlet right-hand-side correction overflowed."
                    );
                }
            } else {
                constrained_entries.push_back(entry);
            }
        }
    }
    for (std::size_t node = 0; node < node_count; ++node) {
        if (prescribed[node].has_value()) {
            rhs[node] = *prescribed[node];
            constrained_entries.push_back({node, node, 1.0});
        }
    }

    const linalg::CsrMatrix matrix = linalg::CsrMatrix::from_triplets(
        node_count, node_count, constrained_entries
    );
    linalg::SolverResult solver = linalg::solve_conjugate_gradient(
        matrix, std::span<const double>(rhs), solver_options
    );
    return {solver.solution, std::move(solver)};
}

} // namespace nova::fem
