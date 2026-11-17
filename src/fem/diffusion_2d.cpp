#include "nova/fem/diffusion_2d.hpp"

#include "nova/fem/triangle_p1.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

struct PhysicalGradient {
    double x{};
    double y{};
};

[[nodiscard]] std::pair<std::size_t, std::size_t> canonical_edge(
    const std::size_t first,
    const std::size_t second
) {
    return std::minmax(first, second);
}

} // namespace

ScalarDiffusion2DResult solve_scalar_diffusion_2d(
    const ScalarDiffusion2DProblem& problem,
    const linalg::SolverOptions& solver_options
) {
    problem.mesh.validate();

    if (!problem.diffusivity) {
        throw std::invalid_argument("A scalar diffusivity function must be provided.");
    }
    if (!problem.source) {
        throw std::invalid_argument("A source function must be provided.");
    }
    if (problem.dirichlet_constraints.empty()) {
        throw std::invalid_argument("At least one Dirichlet constraint is required.");
    }

    const std::size_t node_count = problem.mesh.nodes.size();
    std::vector<std::optional<double>> prescribed_values(node_count);

    for (const DirichletConstraint2D& constraint : problem.dirichlet_constraints) {
        if (constraint.node >= node_count) {
            throw std::out_of_range("Dirichlet constraint references a missing mesh node.");
        }
        if (!std::isfinite(constraint.value)) {
            throw std::invalid_argument("Dirichlet constraint values must be finite.");
        }
        if (prescribed_values[constraint.node].has_value()) {
            throw std::invalid_argument("A mesh node has more than one Dirichlet constraint.");
        }
        prescribed_values[constraint.node] = constraint.value;
    }

    // Build edge incidence so natural boundary data can only be applied
    // to true exterior edges, never to interior edges shared by elements.
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edge_incidence;
    for (const Triangle2D& element : problem.mesh.elements) {
        ++edge_incidence[canonical_edge(element[0], element[1])];
        ++edge_incidence[canonical_edge(element[1], element[2])];
        ++edge_incidence[canonical_edge(element[2], element[0])];
    }

    std::set<std::pair<std::size_t, std::size_t>> prescribed_natural_edges;
    for (const NeumannBoundaryCondition2D& condition : problem.neumann_conditions) {
        const std::size_t first = condition.edge[0];
        const std::size_t second = condition.edge[1];

        if (first >= node_count || second >= node_count) {
            throw std::out_of_range("Neumann edge references a missing mesh node.");
        }
        if (first == second) {
            throw std::invalid_argument("A Neumann edge must have two distinct nodes.");
        }
        if (!condition.outward_flux) {
            throw std::invalid_argument("A Neumann boundary condition requires an outward-flux function.");
        }

        const auto key = canonical_edge(first, second);
        const auto incidence = edge_incidence.find(key);
        if (incidence == edge_incidence.end() || incidence->second != 1U) {
            throw std::invalid_argument("Neumann data must reference an exterior mesh edge.");
        }
        if (!prescribed_natural_edges.insert(key).second) {
            throw std::invalid_argument("A boundary edge has duplicate or conflicting natural conditions.");
        }
    }

    for (const RobinBoundaryCondition2D& condition : problem.robin_conditions) {
        const std::size_t first = condition.edge[0];
        const std::size_t second = condition.edge[1];

        if (first >= node_count || second >= node_count) {
            throw std::out_of_range("Robin edge references a missing mesh node.");
        }
        if (first == second) {
            throw std::invalid_argument("A Robin edge must have two distinct nodes.");
        }
        if (!condition.transfer_coefficient || !condition.ambient_value) {
            throw std::invalid_argument(
                "A Robin condition requires transfer-coefficient and ambient-value functions."
            );
        }

        const auto key = canonical_edge(first, second);
        const auto incidence = edge_incidence.find(key);
        if (incidence == edge_incidence.end() || incidence->second != 1U) {
            throw std::invalid_argument("Robin data must reference an exterior mesh edge.");
        }
        if (!prescribed_natural_edges.insert(key).second) {
            throw std::invalid_argument("A boundary edge has duplicate or conflicting natural conditions.");
        }
    }

    std::vector<linalg::Triplet> matrix_entries;
    const std::size_t maximum_entries = matrix_entries.max_size();
    if (problem.mesh.elements.size() > maximum_entries / 9U) {
        throw std::length_error("The requested 2D finite-element matrix is too large.");
    }
    const std::size_t volume_entry_count = problem.mesh.elements.size() * 9U;
    if (problem.robin_conditions.size() >
        (maximum_entries - volume_entry_count) / 8U) {
        throw std::length_error("The Robin boundary matrix is too large.");
    }
    matrix_entries.reserve(
        volume_entry_count + problem.robin_conditions.size() * 8U
    );

    std::vector<double> right_hand_side(node_count, 0.0);

    for (const Triangle2D& element : problem.mesh.elements) {
        const Point2D& p0 = problem.mesh.nodes[element[0]];
        const Point2D& p1 = problem.mesh.nodes[element[1]];
        const Point2D& p2 = problem.mesh.nodes[element[2]];

        const double determinant =
            (p1.x - p0.x) * (p2.y - p0.y) -
            (p2.x - p0.x) * (p1.y - p0.y);

        const std::array<PhysicalGradient, 3> gradients{{
            {(p1.y - p2.y) / determinant, (p2.x - p1.x) / determinant},
            {(p2.y - p0.y) / determinant, (p0.x - p2.x) / determinant},
            {(p0.y - p1.y) / determinant, (p1.x - p0.x) / determinant}
        }};

        std::array<std::array<double, 3>, 3> local_matrix{};
        std::array<double, 3> local_rhs{};

        for (const TriangleQuadraturePoint& quadrature : triangle_degree_four_quadrature) {
            const std::array<double, 3> shape =
                triangle_p1_shape_values(quadrature.r, quadrature.s);

            const double x =
                shape[0] * p0.x + shape[1] * p1.x + shape[2] * p2.x;
            const double y =
                shape[0] * p0.y + shape[1] * p1.y + shape[2] * p2.y;

            const double coefficient = problem.diffusivity(x, y);
            const double source_value = problem.source(x, y);

            if (!std::isfinite(coefficient) || coefficient <= 0.0) {
                throw std::invalid_argument(
                    "Diffusivity must be finite and strictly positive at every quadrature point."
                );
            }
            if (!std::isfinite(source_value)) {
                throw std::invalid_argument("Source function returned a non-finite value.");
            }

            const double physical_weight =
                std::abs(determinant) * quadrature.weight;

            for (std::size_t i = 0; i < 3; ++i) {
                local_rhs[i] += source_value * shape[i] * physical_weight;

                for (std::size_t j = 0; j < 3; ++j) {
                    const double gradient_product =
                        gradients[i].x * gradients[j].x +
                        gradients[i].y * gradients[j].y;
                    local_matrix[i][j] +=
                        coefficient * gradient_product * physical_weight;
                }
            }
        }

        for (std::size_t i = 0; i < 3; ++i) {
            const std::size_t global_i = element[i];
            right_hand_side[global_i] += local_rhs[i];
            if (!std::isfinite(right_hand_side[global_i])) {
                throw std::overflow_error("Finite-element load-vector accumulation overflowed.");
            }

            for (std::size_t j = 0; j < 3; ++j) {
                matrix_entries.push_back({
                    global_i,
                    element[j],
                    local_matrix[i][j]
                });
            }
        }
    }

    // Two-point Gauss integration on [0,1] for linear boundary edges.
    const double inverse_sqrt_three = 1.0 / std::sqrt(3.0);
    const std::array<double, 2> edge_points{{
        0.5 * (1.0 - inverse_sqrt_three),
        0.5 * (1.0 + inverse_sqrt_three)
    }};

    for (const NeumannBoundaryCondition2D& condition : problem.neumann_conditions) {
        const std::size_t first = condition.edge[0];
        const std::size_t second = condition.edge[1];
        const Point2D& p0 = problem.mesh.nodes[first];
        const Point2D& p1 = problem.mesh.nodes[second];
        const double edge_length = std::hypot(p1.x - p0.x, p1.y - p0.y);

        if (!std::isfinite(edge_length) || edge_length <= 0.0) {
            throw std::invalid_argument("Neumann edge length must be finite and positive.");
        }

        for (const double t : edge_points) {
            const double x = (1.0 - t) * p0.x + t * p1.x;
            const double y = (1.0 - t) * p0.y + t * p1.y;
            const double flux = condition.outward_flux(x, y);

            if (!std::isfinite(flux)) {
                throw std::invalid_argument("Neumann outward-flux function returned a non-finite value.");
            }

            const double integration_weight = 0.5 * edge_length;
            right_hand_side[first] -= flux * (1.0 - t) * integration_weight;
            right_hand_side[second] -= flux * t * integration_weight;

            if (!std::isfinite(right_hand_side[first]) ||
                !std::isfinite(right_hand_side[second])) {
                throw std::overflow_error("Neumann load-vector accumulation overflowed.");
            }
        }
    }

    // Robin convention: q_n = h * (u - u_infinity). Move h*u to the
    // left-hand side and h*u_infinity to the right-hand side.
    for (const RobinBoundaryCondition2D& condition : problem.robin_conditions) {
        const std::size_t first = condition.edge[0];
        const std::size_t second = condition.edge[1];
        const Point2D& p0 = problem.mesh.nodes[first];
        const Point2D& p1 = problem.mesh.nodes[second];
        const double edge_length = std::hypot(p1.x - p0.x, p1.y - p0.y);

        if (!std::isfinite(edge_length) || edge_length <= 0.0) {
            throw std::invalid_argument("Robin edge length must be finite and positive.");
        }

        for (const double t : edge_points) {
            const double x = (1.0 - t) * p0.x + t * p1.x;
            const double y = (1.0 - t) * p0.y + t * p1.y;
            const double transfer = condition.transfer_coefficient(x, y);
            const double ambient = condition.ambient_value(x, y);

            if (!std::isfinite(transfer) || transfer < 0.0) {
                throw std::invalid_argument(
                    "Robin transfer coefficient must be finite and non-negative."
                );
            }
            if (!std::isfinite(ambient)) {
                throw std::invalid_argument("Robin ambient value must be finite.");
            }

            const std::array<double, 2> shape{{1.0 - t, t}};
            const double weight = 0.5 * edge_length;

            right_hand_side[first] += transfer * ambient * shape[0] * weight;
            right_hand_side[second] += transfer * ambient * shape[1] * weight;

            if (!std::isfinite(right_hand_side[first]) ||
                !std::isfinite(right_hand_side[second])) {
                throw std::overflow_error("Robin load-vector accumulation overflowed.");
            }

            matrix_entries.push_back({
                first, first, transfer * shape[0] * shape[0] * weight
            });
            matrix_entries.push_back({
                first, second, transfer * shape[0] * shape[1] * weight
            });
            matrix_entries.push_back({
                second, first, transfer * shape[1] * shape[0] * weight
            });
            matrix_entries.push_back({
                second, second, transfer * shape[1] * shape[1] * weight
            });
        }
    }

    // Strong Dirichlet enforcement: remove constrained columns, correct
    // the free-node RHS, and replace constrained rows with identity rows.
    std::vector<linalg::Triplet> constrained_entries;
    if (node_count > constrained_entries.max_size() - matrix_entries.size()) {
        throw std::length_error("The constrained finite-element matrix is too large.");
    }
    constrained_entries.reserve(matrix_entries.size() + node_count);

    for (const linalg::Triplet& entry : matrix_entries) {
        const bool row_is_constrained = prescribed_values[entry.row].has_value();
        const bool column_is_constrained = prescribed_values[entry.column].has_value();

        if (!row_is_constrained && column_is_constrained) {
            right_hand_side[entry.row] -=
                entry.value * prescribed_values[entry.column].value();
            if (!std::isfinite(right_hand_side[entry.row])) {
                throw std::overflow_error("Dirichlet RHS correction overflowed.");
            }
        } else if (!row_is_constrained && !column_is_constrained) {
            constrained_entries.push_back(entry);
        }
    }

    for (std::size_t node = 0; node < node_count; ++node) {
        if (prescribed_values[node].has_value()) {
            constrained_entries.push_back({node, node, 1.0});
            right_hand_side[node] = prescribed_values[node].value();
        }
    }

    const linalg::CsrMatrix matrix = linalg::CsrMatrix::from_triplets(
        node_count,
        node_count,
        constrained_entries
    );

    linalg::SolverResult linear_solver = linalg::solve_conjugate_gradient(
        matrix,
        std::span<const double>(right_hand_side),
        solver_options
    );

    return {linear_solver.solution, std::move(linear_solver)};
}

} // namespace nova::fem
