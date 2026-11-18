#include "nova/fem/nonlinear_heat_2d.hpp"

#include "nova/fem/triangle_p1.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using EdgeKey = std::pair<std::size_t, std::size_t>;

[[nodiscard]] EdgeKey canonical_edge(
    const std::size_t a, const std::size_t b
) {
    return std::minmax(a, b);
}

struct BoundaryAssignments {
    std::vector<std::optional<double>> dirichlet;
};

[[nodiscard]] BoundaryAssignments validate_boundaries(
    const NonlinearHeat2DProblem& problem,
    const H1FunctionSpace2D& space,
    const double time
) {
    const std::size_t vertex_count = problem.mesh.nodes.size();
    std::map<EdgeKey, std::size_t> incidence;
    for (const Triangle2D& element : problem.mesh.elements) {
        ++incidence[canonical_edge(element[0], element[1])];
        ++incidence[canonical_edge(element[1], element[2])];
        ++incidence[canonical_edge(element[2], element[0])];
    }

    std::set<EdgeKey> assigned;
    const auto check_edge = [&](const BoundaryEdge2D& edge) {
        if (edge[0] >= vertex_count || edge[1] >= vertex_count) {
            throw std::out_of_range("Nonlinear heat edge references a missing mesh vertex.");
        }
        if (edge[0] == edge[1]) {
            throw std::invalid_argument("Nonlinear heat edge must have distinct vertices.");
        }
        const EdgeKey key = canonical_edge(edge[0], edge[1]);
        const auto it = incidence.find(key);
        if (it == incidence.end() || it->second != 1U) {
            throw std::invalid_argument("Nonlinear heat boundary data require an exterior edge.");
        }
        if (!assigned.insert(key).second) {
            throw std::invalid_argument("Overlapping nonlinear heat boundary-edge conditions.");
        }
    };

    BoundaryAssignments result;
    result.dirichlet.resize(space.dimension());
    const auto prescribe = [&](std::size_t index, double value) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Nonlinear heat Dirichlet value must be finite.");
        }
        const auto previous = result.dirichlet[index];
        if (previous.has_value()) {
            const double scale = std::max({
                1.0, std::abs(*previous), std::abs(value)
            });
            if (std::abs(value - *previous) >
                64.0 * std::numeric_limits<double>::epsilon() * scale) {
                throw std::invalid_argument(
                    "Inconsistent nonlinear heat Dirichlet values at a shared DOF."
                );
            }
        }
        result.dirichlet[index] = value;
    };

    for (const auto& bc : problem.dirichlet_edges) {
        check_edge(bc.edge);
        if (!bc.temperature) {
            throw std::invalid_argument("Nonlinear heat Dirichlet function is missing.");
        }
        const Point2D& a = problem.mesh.nodes[bc.edge[0]];
        const Point2D& b = problem.mesh.nodes[bc.edge[1]];
        prescribe(bc.edge[0], bc.temperature(a.x, a.y, time));
        prescribe(bc.edge[1], bc.temperature(b.x, b.y, time));
        if (space.order() == H1TriangleOrder::P2) {
            const std::size_t midpoint =
                space.dofs().edge_dof(bc.edge[0], bc.edge[1]).value();
            prescribe(
                midpoint,
                bc.temperature(
                    0.5 * (a.x + b.x), 0.5 * (a.y + b.y), time
                )
            );
        }
    }
    for (const auto& bc : problem.neumann_edges) {
        check_edge(bc.edge);
        if (!bc.outward_flux) {
            throw std::invalid_argument("Nonlinear heat Neumann function is missing.");
        }
    }
    for (const auto& bc : problem.robin_edges) {
        check_edge(bc.edge);
        if (!bc.transfer_coefficient || !bc.ambient_temperature) {
            throw std::invalid_argument("Nonlinear heat Robin functions are missing.");
        }
    }
    return result;
}

struct EdgeGeometry {
    std::array<std::size_t, 3> dofs{};
    Point2D a{};
    Point2D b{};
    double length{};
};

[[nodiscard]] EdgeGeometry edge_geometry(
    const NonlinearHeat2DProblem& problem,
    const H1FunctionSpace2D& space,
    const BoundaryEdge2D edge
) {
    const Point2D& a = problem.mesh.nodes[edge[0]];
    const Point2D& b = problem.mesh.nodes[edge[1]];
    const double length = std::hypot(b.x - a.x, b.y - a.y);
    if (!std::isfinite(length) || length <= 0.0) {
        throw std::invalid_argument("Nonlinear heat edge length must be positive and finite.");
    }
    const std::size_t midpoint = space.order() == H1TriangleOrder::P2
        ? space.dofs().edge_dof(edge[0], edge[1]).value()
        : 0U;
    return {{edge[0], edge[1], midpoint}, a, b, length};
}

[[nodiscard]] std::array<double, 3> edge_shape(
    const H1TriangleOrder order,
    const double t
) {
    if (order == H1TriangleOrder::P1) {
        return {{1.0 - t, t, 0.0}};
    }
    return {{
        (1.0 - t) * (1.0 - 2.0 * t),
        t * (2.0 * t - 1.0),
        4.0 * t * (1.0 - t)
    }};
}

class InvalidMaterialTrial : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

[[nodiscard]] double valid_positive(const double value, const char* label) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(
            std::string(label) + " must be finite and strictly positive."
        );
    }
    return value;
}

[[nodiscard]] double valid_trial_positive(
    const double value,
    const char* label
) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw InvalidMaterialTrial(
            std::string(label) +
            " is invalid at the current Newton trial temperature."
        );
    }
    return value;
}

} // namespace

NonlinearHeatSystem2D assemble_nonlinear_heat_be_system_2d(
    const NonlinearHeat2DProblem& problem,
    const H1TriangleOrder order,
    const std::span<const double> previous_temperature,
    const std::span<const double> trial_temperature,
    const double time_step,
    const double time
) {
    if (!std::isfinite(time) || !std::isfinite(time_step) ||
        time_step <= 0.0) {
        throw std::invalid_argument("Nonlinear heat time and time step must be finite and positive.");
    }
    if (!problem.volumetric_capacity || !problem.conductivity ||
        !problem.conductivity_temperature_derivative || !problem.source) {
        throw std::invalid_argument(
            "Nonlinear heat capacity, conductivity, derivative and source functions are required."
        );
    }
    const H1FunctionSpace2D space(problem.mesh, order);
    const std::size_t n = space.dimension();
    if (previous_temperature.size() != n || trial_temperature.size() != n) {
        throw std::invalid_argument("Nonlinear heat FEM state dimensions do not match DOFs.");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(previous_temperature[i]) ||
            !std::isfinite(trial_temperature[i])) {
            throw std::invalid_argument("Nonlinear heat state values must be finite.");
        }
    }
    const BoundaryAssignments boundary = validate_boundaries(problem, space, time);
    std::vector<double> residual(n, 0.0);
    std::vector<linalg::Triplet> entries;
    const std::size_t local_size = space.element_dof_count();
    const std::size_t local_entries = local_size * local_size;
    if (problem.mesh.elements.size() > entries.max_size() / local_entries) {
        throw std::length_error("Nonlinear heat Jacobian volume terms are too large.");
    }
    const std::size_t volume_capacity =
        problem.mesh.elements.size() * local_entries;
    if (problem.robin_edges.size() >
        (entries.max_size() - volume_capacity) / 27U) {
        throw std::length_error("Nonlinear heat Robin terms are too large.");
    }
    entries.reserve(volume_capacity + 27U * problem.robin_edges.size());

    for (std::size_t e = 0; e < problem.mesh.elements.size(); ++e) {
        const Triangle2D& cell = problem.mesh.elements[e];
        const auto& ids = space.dofs().element_dofs(e);
        const Point2D& p0 = problem.mesh.nodes[cell[0]];
        const Point2D& p1 = problem.mesh.nodes[cell[1]];
        const Point2D& p2 = problem.mesh.nodes[cell[2]];
        const double j00 = p1.x - p0.x;
        const double j01 = p2.x - p0.x;
        const double j10 = p1.y - p0.y;
        const double j11 = p2.y - p0.y;
        const double det = j00 * j11 - j01 * j10;

        std::array<double, 6> element_residual{};
        std::array<std::array<double, 6>, 6> element_tangent{};

        for (const auto& q : triangle_degree_four_quadrature) {
            const auto basis = evaluate_triangle_h1_basis(order, q.r, q.s);
            const double x = p0.x + j00 * q.r + j01 * q.s;
            const double y = p0.y + j10 * q.r + j11 * q.s;
            std::array<std::array<double, 2>, 6> gradients{};
            double temperature = 0.0;
            double old_temperature = 0.0;
            std::array<double, 2> grad_t{{0.0, 0.0}};
            for (std::size_t i = 0; i < local_size; ++i) {
                const auto& grad = basis.reference_gradients[i];
                gradients[i] = {{
                    (j11 * grad[0] - j10 * grad[1]) / det,
                    (-j01 * grad[0] + j00 * grad[1]) / det
                }};
                temperature += trial_temperature[ids[i]] * basis.values[i];
                old_temperature += previous_temperature[ids[i]] * basis.values[i];
                grad_t[0] += trial_temperature[ids[i]] * gradients[i][0];
                grad_t[1] += trial_temperature[ids[i]] * gradients[i][1];
            }

            const double c = valid_positive(
                problem.volumetric_capacity(x, y), "Volumetric heat capacity"
            );

            const double k = valid_trial_positive(
                problem.conductivity(x, y, temperature),
                "Thermal conductivity"
            );
            const double dk =
                problem.conductivity_temperature_derivative(x, y, temperature);
            if (!std::isfinite(dk)) {
                throw InvalidMaterialTrial(
                    "Thermal conductivity derivative is non-finite at the current Newton trial temperature."
                );
            }

            const double f = problem.source(x, y, time);
            if (!std::isfinite(f)) {
                throw std::invalid_argument("Nonlinear heat source is non-finite.");
            }

            const double weight = std::abs(det) * q.weight;
            for (std::size_t i = 0; i < local_size; ++i) {
                const double grad_dot_test =
                    grad_t[0] * gradients[i][0] +
                    grad_t[1] * gradients[i][1];
                element_residual[i] += weight * (
                    c * (temperature - old_temperature) / time_step *
                        basis.values[i] +
                    k * grad_dot_test - f * basis.values[i]
                );
                for (std::size_t j = 0; j < local_size; ++j) {
                    const double grad_dot_grad =
                        gradients[i][0] * gradients[j][0] +
                        gradients[i][1] * gradients[j][1];
                    element_tangent[i][j] += weight * (
                        c / time_step * basis.values[i] * basis.values[j] +
                        k * grad_dot_grad +
                        dk * basis.values[j] * grad_dot_test
                    );
                }
            }
        }

        for (std::size_t i = 0; i < local_size; ++i) {
            residual[ids[i]] += element_residual[i];
            if (!std::isfinite(residual[ids[i]])) {
                throw std::overflow_error("Nonlinear heat volumetric residual overflowed.");
            }
            for (std::size_t j = 0; j < local_size; ++j) {
                if (!std::isfinite(element_tangent[i][j])) {
                    throw std::overflow_error("Nonlinear heat volumetric tangent overflowed.");
                }
                entries.push_back({
                    ids[i], ids[j], element_tangent[i][j]
                });
            }
        }
    }

    // Three-point Gauss-Legendre integration on the boundary edges.
    constexpr double abscissa = 0.77459666924148337704;
    constexpr std::array<double, 3> nodes{{
        0.5 * (1.0 - abscissa), 0.5, 0.5 * (1.0 + abscissa)
    }};
    constexpr std::array<double, 3> weights{{
        5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0
    }};
    const std::size_t edge_size =
        order == H1TriangleOrder::P2 ? 3U : 2U;

    for (const auto& bc : problem.neumann_edges) {
        const EdgeGeometry edge = edge_geometry(problem, space, bc.edge);
        for (std::size_t q = 0; q < nodes.size(); ++q) {
            const double t = nodes[q];
            const double x = (1.0 - t) * edge.a.x + t * edge.b.x;
            const double y = (1.0 - t) * edge.a.y + t * edge.b.y;
            const double flux = bc.outward_flux(x, y, time);
            if (!std::isfinite(flux)) {
                throw std::invalid_argument("Nonlinear heat Neumann flux must be finite.");
            }
            const auto shape = edge_shape(order, t);
            for (std::size_t i = 0; i < edge_size; ++i) {
                // q_n positive outward; residual = mass + conduction + q_n - source.
                residual[edge.dofs[i]] +=
                    flux * shape[i] * edge.length * weights[q];
            }
        }
    }

    for (const auto& bc : problem.robin_edges) {
        const EdgeGeometry edge = edge_geometry(problem, space, bc.edge);
        for (std::size_t q = 0; q < nodes.size(); ++q) {
            const double t = nodes[q];
            const double x = (1.0 - t) * edge.a.x + t * edge.b.x;
            const double y = (1.0 - t) * edge.a.y + t * edge.b.y;
            const double h = bc.transfer_coefficient(x, y, time);
            const double ambient = bc.ambient_temperature(x, y, time);
            if (!std::isfinite(h) || h < 0.0 || !std::isfinite(ambient)) {
                throw std::invalid_argument(
                    "Nonlinear heat Robin coefficient/ambient must be finite; h nonnegative."
                );
            }
            const auto shape = edge_shape(order, t);
            const double weight = weights[q] * edge.length;
            double edge_temperature = 0.0;
            for (std::size_t j = 0; j < edge_size; ++j) {
                edge_temperature +=
                    shape[j] * trial_temperature[edge.dofs[j]];
            }
            for (std::size_t i = 0; i < edge_size; ++i) {
                residual[edge.dofs[i]] +=
                    h * (edge_temperature - ambient) * shape[i] * weight;
                for (std::size_t j = 0; j < edge_size; ++j) {
                    entries.push_back({
                        edge.dofs[i], edge.dofs[j],
                        h * shape[i] * shape[j] * weight
                    });
                }
            }
        }
    }

    // Strong essential boundary enforcement uses identity rows. Unlike
    // symmetric CG elimination, free rows retain constrained columns
    // so the nonsymmetric Newton Jacobian remains the exact derivative
    // of the residual with respect to the complete trial DOF vector.
    std::vector<linalg::Triplet> constrained_entries;
    if (n > constrained_entries.max_size() - entries.size()) {
        throw std::length_error("Nonlinear heat constrained Jacobian too large.");
    }
    constrained_entries.reserve(entries.size() + n);
    for (const auto& entry : entries) {
        if (!boundary.dirichlet[entry.row].has_value()) {
            constrained_entries.push_back(entry);
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (boundary.dirichlet[i].has_value()) {
            residual[i] = trial_temperature[i] - *boundary.dirichlet[i];
            constrained_entries.push_back({i, i, 1.0});
        }
        if (!std::isfinite(residual[i])) {
            throw std::overflow_error("Nonlinear heat assembled residual is non-finite.");
        }
    }

    return {
        std::move(residual),
        linalg::CsrMatrix::from_triplets(n, n, constrained_entries)
    };
}

NonlinearHeat2DResult solve_nonlinear_heat_2d(
    const NonlinearHeat2DProblem& problem,
    const NonlinearHeat2DOptions& options
) {
    if (!std::isfinite(options.start_time) ||
        !std::isfinite(options.time_step) || options.time_step <= 0.0 ||
        options.steps == 0U) {
        throw std::invalid_argument(
            "Nonlinear heat requires finite start, positive dt and positive step count."
        );
    }
    const double final_time = options.start_time +
        options.time_step * static_cast<double>(options.steps);
    if (!std::isfinite(final_time) || final_time <= options.start_time) {
        throw std::invalid_argument("Nonlinear heat final time must advance and remain finite.");
    }
    if (problem.initial_dof_values.empty() && !problem.initial_temperature) {
        throw std::invalid_argument("Nonlinear heat needs initial DOFs or initial field.");
    }
    const H1FunctionSpace2D space(problem.mesh, options.order);
    std::vector<double> state = problem.initial_dof_values.empty()
        ? space.interpolate(problem.initial_temperature)
        : problem.initial_dof_values;
    if (state.size() != space.dimension()) {
        throw std::invalid_argument("Nonlinear heat initial state DOF mismatch.");
    }
    for (double value : state) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Nonlinear heat initial temperature is non-finite.");
        }
    }
    const auto initial_boundary = validate_boundaries(
        problem, space, options.start_time
    );
    for (std::size_t i = 0; i < state.size(); ++i) {
        if (initial_boundary.dirichlet[i].has_value()) {
            state[i] = *initial_boundary.dirichlet[i];
        }
    }

    NonlinearHeat2DResult result;
    result.final_time = final_time;
    result.newton_steps.reserve(options.steps);

    for (std::size_t step = 1; step <= options.steps; ++step) {
        const double time =
            options.start_time + static_cast<double>(step) * options.time_step;
        if (!std::isfinite(time)) {
            throw std::overflow_error("Nonlinear heat step time overflowed.");
        }
        const auto boundary = validate_boundaries(problem, space, time);
        std::vector<double> starting_guess = state;
        for (std::size_t i = 0; i < starting_guess.size(); ++i) {
            if (boundary.dirichlet[i].has_value()) {
                starting_guess[i] = *boundary.dirichlet[i];
            }
        }
        const auto residual_fn = [&](const std::span<const double> trial) {
            try {
                return assemble_nonlinear_heat_be_system_2d(
                    problem, options.order, state, trial,
                    options.time_step, time
                ).residual;
            } catch (const InvalidMaterialTrial&) {
                // An invalid temperature-dependent material evaluation is a
                // rejected line-search point, not a fatal Newton-step error.
                return std::vector<double>(
                    trial.size(), std::numeric_limits<double>::quiet_NaN()
                );
            }
        };
        const auto jacobian_fn = [&](const std::span<const double> trial) {
            return assemble_nonlinear_heat_be_system_2d(
                problem, options.order, state, trial, options.time_step, time
            ).jacobian;
        };
        nonlinear::NewtonResult newton = nonlinear::solve_newton(
            residual_fn, jacobian_fn,
            std::span<const double>(starting_guess), options.newton_options
        );
        if (!newton.converged) {
            throw std::runtime_error(
                "Nonlinear heat Newton failed at step " +
                std::to_string(step) + ": " + newton.message
            );
        }
        state = newton.solution;
        result.newton_steps.push_back(std::move(newton));
    }

    result.temperature = std::move(state);
    return result;
}

} // namespace nova::fem
