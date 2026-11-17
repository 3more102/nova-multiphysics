#include "nova/fem/poisson_1d.hpp"

#include "nova/linalg/csr_matrix.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {

Poisson1DResult solve_poisson_1d(
    const Poisson1DProblem& problem,
    const linalg::SolverOptions& solver_options
) {
    if (!std::isfinite(problem.length) || problem.length <= 0.0) {
        throw std::invalid_argument("Poisson domain length must be finite and positive.");
    }
    if (problem.elements == 0) {
        throw std::invalid_argument("Poisson mesh must contain at least one element.");
    }
    if (!std::isfinite(problem.conductivity) || problem.conductivity <= 0.0) {
        throw std::invalid_argument("Poisson conductivity must be finite and positive.");
    }
    if (!std::isfinite(problem.left_value) || !std::isfinite(problem.right_value)) {
        throw std::invalid_argument("Poisson Dirichlet values must be finite.");
    }
    if (!problem.source) {
        throw std::invalid_argument("Poisson source function must be provided.");
    }

    const std::size_t node_count = problem.elements + 1;
    const std::size_t last_node = node_count - 1;
    const double element_length = problem.length / static_cast<double>(problem.elements);

    std::vector<double> coordinates(node_count, 0.0);
    for (std::size_t node = 0; node < node_count; ++node) {
        coordinates[node] = static_cast<double>(node) * element_length;
    }
    coordinates.back() = problem.length;

    std::vector<linalg::Triplet> triplets;
    triplets.reserve(4 * problem.elements);
    std::vector<double> right_hand_side(node_count, 0.0);

    const double inverse_sqrt_three = 1.0 / std::sqrt(3.0);
    const std::array<double, 2> quadrature_points{
        0.5 * (1.0 - inverse_sqrt_three),
        0.5 * (1.0 + inverse_sqrt_three)
    };

    for (std::size_t element = 0; element < problem.elements; ++element) {
        const std::size_t node_0 = element;
        const std::size_t node_1 = element + 1;
        const double x_0 = coordinates[node_0];

        const double stiffness = problem.conductivity / element_length;
        triplets.push_back({node_0, node_0, stiffness});
        triplets.push_back({node_0, node_1, -stiffness});
        triplets.push_back({node_1, node_0, -stiffness});
        triplets.push_back({node_1, node_1, stiffness});

        for (const double xi : quadrature_points) {
            const double physical_x = x_0 + xi * element_length;
            const double source_value = problem.source(physical_x);
            if (!std::isfinite(source_value)) {
                throw std::invalid_argument("Poisson source returned a non-finite value.");
            }

            const double shape_0 = 1.0 - xi;
            const double shape_1 = xi;
            const double integration_weight = 0.5 * element_length;

            right_hand_side[node_0] += source_value * shape_0 * integration_weight;
            right_hand_side[node_1] += source_value * shape_1 * integration_weight;
        }
    }

    // Strong Dirichlet enforcement by column elimination and RHS correction.
    std::vector<linalg::Triplet> constrained_triplets;
    constrained_triplets.reserve(triplets.size() + 2);

    for (const linalg::Triplet& entry : triplets) {
        const bool row_is_boundary = entry.row == 0 || entry.row == last_node;
        const bool column_is_boundary = entry.column == 0 || entry.column == last_node;

        if (!row_is_boundary && column_is_boundary) {
            const double prescribed_value =
                entry.column == 0 ? problem.left_value : problem.right_value;
            right_hand_side[entry.row] -= entry.value * prescribed_value;
        } else if (!row_is_boundary && !column_is_boundary) {
            constrained_triplets.push_back(entry);
        }
    }

    constrained_triplets.push_back({0, 0, 1.0});
    constrained_triplets.push_back({last_node, last_node, 1.0});
    right_hand_side[0] = problem.left_value;
    right_hand_side[last_node] = problem.right_value;

    const linalg::CsrMatrix matrix = linalg::CsrMatrix::from_triplets(
        node_count,
        node_count,
        constrained_triplets
    );

    linalg::SolverResult linear_solver = linalg::solve_conjugate_gradient(
        matrix,
        std::span<const double>(right_hand_side),
        solver_options
    );

    return {
        std::move(coordinates),
        linear_solver.solution,
        std::move(linear_solver)
    };
}

} // namespace nova::fem
