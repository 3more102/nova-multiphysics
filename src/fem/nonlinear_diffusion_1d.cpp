#include "nova/fem/nonlinear_diffusion_1d.hpp"

#include "nova/linalg/csr_matrix.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {

NonlinearDiffusion1DResult solve_nonlinear_diffusion_1d(
    const NonlinearDiffusion1DProblem& problem,
    const std::span<const double> initial_guess,
    const nonlinear::NewtonOptions& newton_options
) {
    if (!std::isfinite(problem.length) || problem.length <= 0.0) {
        throw std::invalid_argument("Nonlinear diffusion domain length must be finite and positive.");
    }
    if (problem.elements == 0 ||
        problem.elements == std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument("Nonlinear diffusion mesh must contain a valid positive element count.");
    }
    if (!problem.source) {
        throw std::invalid_argument("Nonlinear diffusion source function must be provided.");
    }
    if (!std::isfinite(problem.left_value) || !std::isfinite(problem.right_value)) {
        throw std::invalid_argument("Nonlinear diffusion Dirichlet values must be finite.");
    }

    const std::size_t node_count = problem.elements + 1;
    if (!initial_guess.empty() && initial_guess.size() != node_count) {
        throw std::invalid_argument("Nonlinear diffusion initial-guess size must equal the mesh node count.");
    }

    const double element_length =
        problem.length / static_cast<double>(problem.elements);
    if (!std::isfinite(element_length) || element_length <= 0.0) {
        throw std::invalid_argument("Nonlinear diffusion element length is outside the supported numeric range.");
    }

    std::vector<double> coordinates(node_count, 0.0);
    for (std::size_t node = 0; node < node_count; ++node) {
        coordinates[node] = static_cast<double>(node) * element_length;
    }
    coordinates.back() = problem.length;

    std::vector<double> starting_values(node_count, 0.0);
    if (!initial_guess.empty()) {
        starting_values.assign(initial_guess.begin(), initial_guess.end());
    }
    for (const double value : starting_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Nonlinear diffusion initial guess must contain only finite values.");
        }
    }

    // Keep the constrained degrees of freedom fixed throughout Newton updates.
    starting_values.front() = problem.left_value;
    starting_values.back() = problem.right_value;

    constexpr double root_three_fifths = 0.77459666924148337704;
    const std::array<double, 3> quadrature_points{{
        0.5 * (1.0 - root_three_fifths),
        0.5,
        0.5 * (1.0 + root_three_fifths)
    }};
    const std::array<double, 3> quadrature_weights{{
        5.0 / 18.0,
        4.0 / 9.0,
        5.0 / 18.0
    }};

    const nonlinear::NonlinearResidualFunction residual_function =
        [&](const std::span<const double> values) {
            std::vector<double> residual(node_count, 0.0);

            for (std::size_t element = 0; element < problem.elements; ++element) {
                const std::size_t node_0 = element;
                const std::size_t node_1 = element + 1;
                const double x_0 = coordinates[node_0];
                const double h = coordinates[node_1] - x_0;
                const double field_gradient =
                    (values[node_1] - values[node_0]) / h;
                const std::array<double, 2> shape_gradients{{
                    -1.0 / h,
                    1.0 / h
                }};

                for (std::size_t q = 0; q < quadrature_points.size(); ++q) {
                    const double xi = quadrature_points[q];
                    const std::array<double, 2> shape{{1.0 - xi, xi}};
                    const double x = x_0 + xi * h;
                    const double field_value =
                        shape[0] * values[node_0] +
                        shape[1] * values[node_1];
                    const double conductivity =
                        1.0 + field_value * field_value;
                    const double source_value = problem.source(x);

                    if (!std::isfinite(conductivity) ||
                        !std::isfinite(source_value)) {
                        throw std::invalid_argument(
                            "Nonlinear diffusion material or source evaluation became non-finite."
                        );
                    }

                    const double weight = h * quadrature_weights[q];
                    for (std::size_t i = 0; i < 2; ++i) {
                        residual[element + i] +=
                            (conductivity * field_gradient * shape_gradients[i] -
                             source_value * shape[i]) * weight;
                    }
                }
            }

            residual.front() = values.front() - problem.left_value;
            residual.back() = values.back() - problem.right_value;
            return residual;
        };

    const nonlinear::NonlinearJacobianFunction jacobian_function =
        [&](const std::span<const double> values) {
            std::vector<linalg::Triplet> entries;
            entries.reserve(4 * problem.elements + 2);
            const std::size_t last_node = node_count - 1;

            for (std::size_t element = 0; element < problem.elements; ++element) {
                const std::size_t node_0 = element;
                const std::size_t node_1 = element + 1;
                const double x_0 = coordinates[node_0];
                const double h = coordinates[node_1] - x_0;
                const double field_gradient =
                    (values[node_1] - values[node_0]) / h;
                const std::array<double, 2> shape_gradients{{
                    -1.0 / h,
                    1.0 / h
                }};

                for (std::size_t q = 0; q < quadrature_points.size(); ++q) {
                    const double xi = quadrature_points[q];
                    const std::array<double, 2> shape{{1.0 - xi, xi}};
                    const double field_value =
                        shape[0] * values[node_0] +
                        shape[1] * values[node_1];
                    const double conductivity =
                        1.0 + field_value * field_value;
                    const double conductivity_derivative = 2.0 * field_value;
                    const double weight = h * quadrature_weights[q];

                    if (!std::isfinite(conductivity) ||
                        !std::isfinite(conductivity_derivative)) {
                        throw std::invalid_argument(
                            "Nonlinear diffusion Jacobian evaluation became non-finite."
                        );
                    }

                    for (std::size_t i = 0; i < 2; ++i) {
                        const std::size_t global_row = element + i;

                        // Strong Dirichlet rows are replaced by identity rows below.
                        if (global_row == 0 || global_row == last_node) {
                            continue;
                        }

                        for (std::size_t j = 0; j < 2; ++j) {
                            const double tangent =
                                (conductivity_derivative * shape[j] *
                                     field_gradient * shape_gradients[i] +
                                 conductivity * shape_gradients[i] *
                                     shape_gradients[j]) * weight;

                            entries.push_back({
                                global_row,
                                element + j,
                                tangent
                            });
                        }
                    }
                }
            }

            entries.push_back({0, 0, 1.0});
            entries.push_back({last_node, last_node, 1.0});

            return linalg::CsrMatrix::from_triplets(
                node_count,
                node_count,
                entries
            );
        };

    nonlinear::NewtonResult nonlinear_solver = nonlinear::solve_newton(
        residual_function,
        jacobian_function,
        std::span<const double>(starting_values),
        newton_options
    );

    return {
        std::move(coordinates),
        std::move(nonlinear_solver)
    };
}

} // namespace nova::fem
