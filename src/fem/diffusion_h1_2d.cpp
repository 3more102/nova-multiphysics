#include "nova/fem/diffusion_2d.hpp"

#include "nova/fem/h1_space_2d.hpp"
#include "nova/fem/triangle_p1.hpp"
#include "diffusion_nullspace.hpp"
#include "constrained_diffusion.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using Edge = std::pair<std::size_t, std::size_t>;

[[nodiscard]] Edge canonical_edge(std::size_t a, std::size_t b) {
    return std::minmax(a, b);
}

struct EdgeQuadraturePoint {
    double t{};
    double weight{};
};

[[nodiscard]] std::array<EdgeQuadraturePoint, 3> edge_quadrature() {
    // Gauss-Legendre degree five on [0,1]; needed for P2 Robin products.
    const double offset = 0.5 * std::sqrt(3.0 / 5.0);
    return {{{0.5 - offset, 5.0 / 18.0},
             {0.5,          8.0 / 18.0},
             {0.5 + offset, 5.0 / 18.0}}};
}

[[nodiscard]] std::array<double, 3> edge_shape(
    H1TriangleOrder order, double t
) {
    if (order == H1TriangleOrder::P1) {
        return {1.0 - t, t, 0.0};
    }
    // Quadratic Lagrange edge: first endpoint, second endpoint, midpoint.
    return {(1.0 - t) * (1.0 - 2.0 * t),
            t * (2.0 * t - 1.0),
            4.0 * t * (1.0 - t)};
}

void prescribe(
    std::vector<std::optional<double>>& prescribed,
    std::size_t dof, double value
) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("Dirichlet values must be finite.");
    }
    if (prescribed[dof].has_value()) {
        const double old = *prescribed[dof];
        const double scale = std::max({1.0, std::abs(old), std::abs(value)});
        if (std::abs(old - value) > 64.0 * std::numeric_limits<double>::epsilon() * scale) {
            throw std::invalid_argument("Conflicting Dirichlet values at a shared DOF.");
        }
    }
    prescribed[dof] = value;
}

} // namespace

ScalarDiffusion2DResult solve_scalar_diffusion_h1_2d(
    const ScalarDiffusion2DProblem& problem,
    H1TriangleOrder order,
    const linalg::SolverOptions& solver_options
) {
    problem.mesh.validate();
    const H1FunctionSpace2D space(problem.mesh, order);
    const H1DofManager2D& dofs = space.dofs();
    const std::size_t vertex_count = problem.mesh.nodes.size();
    const std::size_t n = space.dimension();
    const std::size_t local_n = space.element_dof_count();
    detail::DiffusionNullspace nullspace(n, problem.mesh.elements.size(),
        [&](std::size_t e) {
            const auto& local = dofs.element_dofs(e);
            return std::span<const std::size_t>(local.dofs.data(), local.count);
        });

    if (!problem.diffusivity || !problem.source) {
        throw std::invalid_argument("Diffusivity and source functions are required.");
    }
     // Explicit nodal constraints retain the original P1 API. For P2, use
    // Dirichlet edges to also constrain their shared midpoint DOFs.
    std::vector<std::optional<double>> prescribed(n);
    std::set<std::size_t> explicit_nodes;
    for (const DirichletConstraint2D& constraint : problem.dirichlet_constraints) {
        if (constraint.node >= vertex_count) {
            throw std::out_of_range("Dirichlet node is not a mesh vertex.");
        }
        if (!explicit_nodes.insert(constraint.node).second) {
            throw std::invalid_argument("Duplicate explicit Dirichlet node.");
        }
        prescribe(prescribed, constraint.node, constraint.value);
    }

    std::map<Edge, std::size_t> edge_incidence;
    for (const Triangle2D& cell : problem.mesh.elements) {
        ++edge_incidence[canonical_edge(cell[0], cell[1])];
        ++edge_incidence[canonical_edge(cell[1], cell[2])];
        ++edge_incidence[canonical_edge(cell[2], cell[0])];
    }
    const auto validate_exterior = [&](const BoundaryEdge2D& edge) {
        if (edge[0] >= vertex_count || edge[1] >= vertex_count) {
            throw std::out_of_range("Boundary edge references a missing vertex.");
        }
        if (edge[0] == edge[1]) {
            throw std::invalid_argument("Boundary edge endpoints must be distinct.");
        }
        const Edge key = canonical_edge(edge[0], edge[1]);
        const auto it = edge_incidence.find(key);
        if (it == edge_incidence.end() || it->second != 1U) {
            throw std::invalid_argument("Boundary conditions require an exterior mesh edge.");
        }
        return key;
    };

    std::set<Edge> dirichlet_edges;
    for (const DirichletBoundaryCondition2D& condition : problem.dirichlet_edges) {
        const Edge key = validate_exterior(condition.edge);
        if (!condition.value) {
            throw std::invalid_argument("Dirichlet edge needs a value function.");
        }
        if (!dirichlet_edges.insert(key).second) {
            throw std::invalid_argument("Duplicate Dirichlet boundary edge.");
        }

        const Point2D& a = problem.mesh.nodes[condition.edge[0]];
        const Point2D& b = problem.mesh.nodes[condition.edge[1]];
        prescribe(prescribed, condition.edge[0], condition.value(a.x, a.y));
        prescribe(prescribed, condition.edge[1], condition.value(b.x, b.y));
        if (order == H1TriangleOrder::P2) {
            prescribe(
                prescribed, dofs.edge_dof(condition.edge[0], condition.edge[1]).value(),
                condition.value(0.5 * (a.x + b.x), 0.5 * (a.y + b.y))
            );
        }
    }

    std::set<Edge> natural_edges;
    for (const NeumannBoundaryCondition2D& condition : problem.neumann_conditions) {
        const Edge key = validate_exterior(condition.edge);
        if (!condition.outward_flux) {
            throw std::invalid_argument("Neumann outward flux function is required.");
        }
        if (dirichlet_edges.count(key) != 0U || !natural_edges.insert(key).second) {
            throw std::invalid_argument("Duplicate or incompatible boundary-edge conditions.");
        }
    }
    for (const RobinBoundaryCondition2D& condition : problem.robin_conditions) {
        const Edge key = validate_exterior(condition.edge);
        if (!condition.transfer_coefficient || !condition.ambient_value) {
            throw std::invalid_argument("Robin coefficient and ambient functions are required.");
        }
        if (dirichlet_edges.count(key) != 0U || !natural_edges.insert(key).second) {
            throw std::invalid_argument("Duplicate or incompatible boundary-edge conditions.");
        }
    }

    std::vector<linalg::Triplet> entries;
    const std::size_t volume_per_element = local_n * local_n;
    const std::size_t edge_per_robin = 3U * local_n * local_n;
    const std::size_t capacity = entries.max_size();
    if (problem.mesh.elements.size() > capacity / volume_per_element) {
        throw std::length_error("Finite-element volume matrix is too large.");
    }
    const std::size_t volume_count = problem.mesh.elements.size() * volume_per_element;
    if (problem.robin_conditions.size() > (capacity - volume_count) / edge_per_robin) {
        throw std::length_error("Finite-element Robin matrix is too large.");
    }
    entries.reserve(volume_count + problem.robin_conditions.size() * edge_per_robin);
    std::vector<double> rhs(n, 0.0);

    for (std::size_t e = 0; e < problem.mesh.elements.size(); ++e) {
        const Triangle2D& cell = problem.mesh.elements[e];
        const auto& local_dofs = dofs.element_dofs(e);
        const Point2D& p0 = problem.mesh.nodes[cell[0]];
        const Point2D& p1 = problem.mesh.nodes[cell[1]];
        const Point2D& p2 = problem.mesh.nodes[cell[2]];
        const double j00 = p1.x - p0.x;
        const double j01 = p2.x - p0.x;
        const double j10 = p1.y - p0.y;
        const double j11 = p2.y - p0.y;
        const double det = j00 * j11 - j01 * j10;
        // Exact integrals of the P1/P2 basis functions on this triangle.
        // P2 vertex basis integrals are zero; each edge basis integrates to A/3.
        const long double one_third_area = std::abs(static_cast<long double>(det)) / 6.0L;
        const std::size_t first_weighted = order == H1TriangleOrder::P2 ? 3U : 0U;
        for (std::size_t i = first_weighted; i < local_n; ++i) {
            nullspace.add_basis_integral(local_dofs[i], one_third_area);
        }
        std::array<std::array<double, 6>, 6> local_matrix{};
        std::array<double, 6> local_rhs{};

        for (const TriangleQuadraturePoint& quad : triangle_degree_four_quadrature) {
            const TriangleH1BasisEvaluation basis =
                evaluate_triangle_h1_basis(order, quad.r, quad.s);
            const double x = p0.x + j00 * quad.r + j01 * quad.s;
            const double y = p0.y + j10 * quad.r + j11 * quad.s;
            const double k = problem.diffusivity(x, y);
            const double f = problem.source(x, y);
            if (!std::isfinite(k) || k <= 0.0) {
                throw std::invalid_argument("Diffusivity must be finite and positive.");
            }
            if (!std::isfinite(f)) {
                throw std::invalid_argument("Source must be finite.");
            }
            const double weight = std::abs(det) * quad.weight;
            std::array<std::array<double, 2>, 6> physical_gradients{};
            for (std::size_t i = 0; i < local_n; ++i) {
                const double dr = basis.reference_gradients[i][0];
                const double ds = basis.reference_gradients[i][1];
                physical_gradients[i] = {
                    (j11 * dr - j10 * ds) / det,
                    (-j01 * dr + j00 * ds) / det
                };
            }
            for (std::size_t i = 0; i < local_n; ++i) {
                const double load = f * basis.values[i] * weight;
                local_rhs[i] += load;
                nullspace.record_load(local_dofs[i], load);
                for (std::size_t j = 0; j < local_n; ++j) {
                    const double dot =
                        physical_gradients[i][0] * physical_gradients[j][0] +
                        physical_gradients[i][1] * physical_gradients[j][1];
                    local_matrix[i][j] += k * dot * weight;
                }
            }
        }

        for (std::size_t i = 0; i < local_n; ++i) {
            const std::size_t global_i = local_dofs[i];
            rhs[global_i] += local_rhs[i];
            if (!std::isfinite(rhs[global_i])) {
                throw std::overflow_error("Element load-vector accumulation overflowed.");
            }
            for (std::size_t j = 0; j < local_n; ++j) {
                entries.push_back({global_i, local_dofs[j], local_matrix[i][j]});
            }
        }
    }

    const auto quadrature = edge_quadrature();
    const auto edge_data = [&](const BoundaryEdge2D& edge) {
        const Point2D& a = problem.mesh.nodes[edge[0]];
        const Point2D& b = problem.mesh.nodes[edge[1]];
        const double length = std::hypot(b.x - a.x, b.y - a.y);
        if (!std::isfinite(length) || length <= 0.0) {
            throw std::invalid_argument("Boundary edge length must be finite and positive.");
        }
        const std::array<std::size_t, 3> ids{{
            edge[0], edge[1],
            order == H1TriangleOrder::P2
                ? dofs.edge_dof(edge[0], edge[1]).value() : 0U
        }};
        return std::pair{ids, length};
    };
    const std::size_t edge_n = order == H1TriangleOrder::P2 ? 3U : 2U;

    // PDE convention: -div(k grad u) = f; outward flux q=-k grad u.n.
    // Weak Neumann contribution is -integral(q * v) on the RHS.
    for (const NeumannBoundaryCondition2D& condition : problem.neumann_conditions) {
        const auto [ids, length] = edge_data(condition.edge);
        const Point2D& a = problem.mesh.nodes[condition.edge[0]];
        const Point2D& b = problem.mesh.nodes[condition.edge[1]];
        for (const EdgeQuadraturePoint& q : quadrature) {
            const double x = (1.0 - q.t) * a.x + q.t * b.x;
            const double y = (1.0 - q.t) * a.y + q.t * b.y;
            const double flux = condition.outward_flux(x, y);
            if (!std::isfinite(flux)) {
                throw std::invalid_argument("Neumann flux must be finite.");
            }
            const auto phi = edge_shape(order, q.t);
            for (std::size_t i = 0; i < edge_n; ++i) {
                const double load = flux * phi[i] * length * q.weight;
                rhs[ids[i]] -= load;
                nullspace.record_load(ids[i], load);
                if (!std::isfinite(rhs[ids[i]])) {
                    throw std::overflow_error("Neumann load-vector accumulation overflowed.");
                }
            }
        }
    }

    // q_n = h * (u - u_infinity): +integral(h u v) on LHS.
    for (const RobinBoundaryCondition2D& condition : problem.robin_conditions) {
        const auto [ids, length] = edge_data(condition.edge);
        const Point2D& a = problem.mesh.nodes[condition.edge[0]];
        const Point2D& b = problem.mesh.nodes[condition.edge[1]];
        for (const EdgeQuadraturePoint& q : quadrature) {
            const double x = (1.0 - q.t) * a.x + q.t * b.x;
            const double y = (1.0 - q.t) * a.y + q.t * b.y;
            const double h = condition.transfer_coefficient(x, y);
            const double ambient = condition.ambient_value(x, y);
            if (!std::isfinite(h) || h < 0.0) {
                throw std::invalid_argument("Robin transfer coefficient must be nonnegative and finite.");
            }
            if (!std::isfinite(ambient)) {
                throw std::invalid_argument("Robin ambient value must be finite.");
            }
            const auto phi = edge_shape(order, q.t);
            if (h > 0.0) {
                nullspace.anchor(ids[0]);
            }
            const double weight = length * q.weight;
            for (std::size_t i = 0; i < edge_n; ++i) {
                rhs[ids[i]] += h * ambient * phi[i] * weight;
                if (!std::isfinite(rhs[ids[i]])) {
                    throw std::overflow_error("Robin load-vector accumulation overflowed.");
                }
                for (std::size_t j = 0; j < edge_n; ++j) {
                    entries.push_back({ids[i], ids[j], h * phi[i] * phi[j] * weight});
                }
            }
        }
    }

    nullspace.impose_gauges(rhs, prescribed);

    auto solver = detail::solve_constrained_diffusion(
        std::move(entries), std::move(rhs), prescribed, solver_options);
    return {nullspace.physical_values(solver.solution), std::move(solver)};
}

} // namespace nova::fem
