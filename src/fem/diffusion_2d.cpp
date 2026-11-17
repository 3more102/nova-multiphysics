#include "nova/fem/diffusion_2d.hpp"

#include "nova/fem/triangle_p1.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

struct PhysicalGradient {
    double x{};
    double y{};
};

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

    std::vector<linalg::Triplet> matrix_entries;
    matrix_entries.reserve(problem.mesh.elements.size() * 9);
    std::vector<double> right_hand_side(node_count, 0.0);

    for (const Triangle2D& element : problem.mesh.elements) {
        const Point2D& p0 = problem.mesh.nodes[element[0]];
        const Point2D& p1 = problem.mesh.nodes[element[1]];
        const Point2D& p2 = problem.mesh.nodes[element[2]];

        const double determinant =
            (p1.x - p0.x) * (p2.y - p0.y) -
            (p2.x - p0.x) * (p1.y - p0.y);
        const double area = 0.5 * std::abs(determinant);

        const std::array<PhysicalGradient, 3> gradients{{
            {(p1.y - p2.y) / determinant, (p2.x - p1.x) / determinant},
            {(p2.y - p0.y) / determinant, (p0.x - p2.x) / determinant},
            {(p0.y - p1.y) / determinant, (p1.x - p0.x) / determinant}
        }};

        std::array<std::array<double, 3>, 3> local_matrix{};
        std::array<double, 3> local_rhs{};

        for (const TriangleQuadraturePoint& quadrature : triangle_degree_two_quadrature) {
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

            const double physical_weight = std::abs(determinant) * quadrature.weight;

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

            for (std::size_t j = 0; j < 3; ++j) {
                matrix_entries.push_back({
                    global_i,
                    element[j],
                    local_matrix[i][j]
                });
            }
        }
    }

    // Strong Dirichlet enforcement: remove constrained columns, correct
    // the free-node RHS, and replace constrained rows with identity rows.
    std::vector<linalg::Triplet> constrained_entries;
    constrained_entries.reserve(matrix_entries.size() + node_count);

    for (const linalg::Triplet& entry : matrix_entries) {
        const bool row_is_constrained = prescribed_values[entry.row].has_value();
        const bool column_is_constrained = prescribed_values[entry.column].has_value();

        if (!row_is_constrained && column_is_constrained) {
            right_hand_side[entry.row] -=
                entry.value * prescribed_values[entry.column].value();
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
