#include "nova/fem/transient_heat_2d.hpp"

#include "nova/fem/triangle_p1.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using Edge = std::pair<std::size_t, std::size_t>;

[[nodiscard]] Edge edge_key(std::size_t a, std::size_t b) {
    return std::minmax(a, b);
}

[[nodiscard]] double checked_positive(double value, const char* name) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be finite and positive.");
    }
    return value;
}

struct ElementMapping {
    double j00{}, j01{}, j10{}, j11{}, determinant{};
    Point2D origin{};
};

[[nodiscard]] ElementMapping map_element(
    const TriangleMesh2D& mesh,
    const Triangle2D& element
) {
    const Point2D& a = mesh.nodes[element[0]];
    const Point2D& b = mesh.nodes[element[1]];
    const Point2D& c = mesh.nodes[element[2]];
    const double j00 = b.x - a.x;
    const double j01 = c.x - a.x;
    const double j10 = b.y - a.y;
    const double j11 = c.y - a.y;
    return {j00, j01, j10, j11, j00 * j11 - j01 * j10, a};
}

[[nodiscard]] Point2D physical_point(
    const ElementMapping& map, double r, double s
) {
    return {
        map.origin.x + map.j00 * r + map.j01 * s,
        map.origin.y + map.j10 * r + map.j11 * s
    };
}

[[nodiscard]] std::vector<double> assemble_load(
    const H1FunctionSpace2D& space,
    const TransientScalarField2D& source,
    double time
) {
    std::vector<double> load(space.dimension(), 0.0);
    const TriangleMesh2D& mesh = space.mesh();
    const std::size_t count = space.element_dof_count();
    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        const ElementMapping map = map_element(mesh, mesh.elements[e]);
        const auto& ids = space.dofs().element_dofs(e);
        for (const TriangleQuadraturePoint& q : triangle_degree_four_quadrature) {
            const Point2D point = physical_point(map, q.r, q.s);
            const double f = source(point.x, point.y, time);
            if (!std::isfinite(f)) {
                throw std::invalid_argument("Transient source must return finite values.");
            }
            const auto basis = evaluate_triangle_h1_basis(space.order(), q.r, q.s);
            const double weight = std::abs(map.determinant) * q.weight;
            for (std::size_t i = 0; i < count; ++i) {
                load[ids[i]] += f * basis.values[i] * weight;
                if (!std::isfinite(load[ids[i]])) {
                    throw std::overflow_error("Transient load-vector accumulation overflowed.");
                }
            }
        }
    }
    return load;
}

struct BoundaryData {
    std::vector<Point2D> coordinates;
    std::vector<std::vector<TransientScalarField2D>> functions;
    std::vector<bool> constrained;
};

[[nodiscard]] BoundaryData prepare_boundary_data(
    const TransientHeat2DProblem& problem,
    const H1FunctionSpace2D& space
) {
    const TriangleMesh2D& mesh = problem.mesh;
    const auto& dofs = space.dofs();
    const std::size_t vertex_count = mesh.nodes.size();
    BoundaryData data;
    data.coordinates = mesh.nodes;
    if (space.order() == H1TriangleOrder::P2) {
        for (const Edge& edge : dofs.edge_keys()) {
            const Point2D& a = mesh.nodes[edge.first];
            const Point2D& b = mesh.nodes[edge.second];
            data.coordinates.push_back({
                0.5 * (a.x + b.x), 0.5 * (a.y + b.y)
            });
        }
    }
    data.functions.resize(space.dimension());
    data.constrained.assign(space.dimension(), false);

    std::map<Edge, std::size_t> incidence;
    for (const Triangle2D& cell : mesh.elements) {
        ++incidence[edge_key(cell[0], cell[1])];
        ++incidence[edge_key(cell[1], cell[2])];
        ++incidence[edge_key(cell[2], cell[0])];
    }
    std::set<Edge> used;
    for (const TransientDirichletEdge2D& condition : problem.dirichlet_edges) {
        const auto [first, second] = condition.edge;
        if (first >= vertex_count || second >= vertex_count) {
            throw std::out_of_range("Transient Dirichlet edge references a missing node.");
        }
        if (first == second) {
            throw std::invalid_argument("Transient Dirichlet edge has identical endpoints.");
        }
        const Edge edge = edge_key(first, second);
        const auto it = incidence.find(edge);
        if (it == incidence.end() || it->second != 1U) {
            throw std::invalid_argument("Transient Dirichlet edge must be exterior.");
        }
        if (!condition.temperature) {
            throw std::invalid_argument("Transient Dirichlet temperature function is required.");
        }
        if (!used.insert(edge).second) {
            throw std::invalid_argument("Duplicate transient Dirichlet edge.");
        }
        data.functions[first].push_back(condition.temperature);
        data.functions[second].push_back(condition.temperature);
        if (space.order() == H1TriangleOrder::P2) {
            const std::size_t midpoint_dof = dofs.edge_dof(first, second).value();
            data.functions[midpoint_dof].push_back(condition.temperature);
        }
    }
    for (std::size_t i = 0; i < data.functions.size(); ++i) {
        data.constrained[i] = !data.functions[i].empty();
    }
    return data;
}

[[nodiscard]] std::vector<double> boundary_values(
    const BoundaryData& boundary, double time
) {
    std::vector<double> result(boundary.functions.size(), 0.0);
    for (std::size_t i = 0; i < boundary.functions.size(); ++i) {
        if (!boundary.constrained[i]) {
            continue;
        }
        const Point2D& p = boundary.coordinates[i];
        const double first = boundary.functions[i][0](p.x, p.y, time);
        if (!std::isfinite(first)) {
            throw std::invalid_argument("Transient Dirichlet value must be finite.");
        }
        for (std::size_t j = 1; j < boundary.functions[i].size(); ++j) {
            const double other = boundary.functions[i][j](p.x, p.y, time);
            if (!std::isfinite(other)) {
                throw std::invalid_argument("Transient Dirichlet value must be finite.");
            }
            const double scale = std::max({1.0, std::abs(first), std::abs(other)});
            if (std::abs(first - other) >
                64.0 * std::numeric_limits<double>::epsilon() * scale) {
                throw std::invalid_argument("Conflicting transient Dirichlet values.");
            }
        }
        result[i] = first;
    }
    return result;
}

} // namespace

ScalarH1Operators2D assemble_scalar_h1_operators_2d(
    const TriangleMesh2D& mesh,
    const H1TriangleOrder order,
    const SpatialScalarField2D& volumetric_capacity,
    const SpatialScalarField2D& conductivity
) {
    if (!volumetric_capacity || !conductivity) {
        throw std::invalid_argument("Volumetric capacity and conductivity functions are required.");
    }
    const H1FunctionSpace2D space(mesh, order);
    const std::size_t n = space.dimension();
    const std::size_t count = space.element_dof_count();
    const std::size_t count_squared = count * count;
    std::vector<linalg::Triplet> mass_entries;
    std::vector<linalg::Triplet> stiffness_entries;
    if (mesh.elements.size() > mass_entries.max_size() / count_squared ||
        mesh.elements.size() > stiffness_entries.max_size() / count_squared) {
        throw std::length_error("Scalar H1 operator matrix is too large.");
    }
    const std::size_t reserve_count = mesh.elements.size() * count_squared;
    mass_entries.reserve(reserve_count);
    stiffness_entries.reserve(reserve_count);

    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        const ElementMapping map = map_element(mesh, mesh.elements[e]);
        const auto& ids = space.dofs().element_dofs(e);
        std::array<std::array<double, 6>, 6> mass{};
        std::array<std::array<double, 6>, 6> stiffness{};

        for (const TriangleQuadraturePoint& q : triangle_degree_four_quadrature) {
            const Point2D point = physical_point(map, q.r, q.s);
            const double c = checked_positive(
                volumetric_capacity(point.x, point.y), "Volumetric capacity"
            );
            const double k = checked_positive(
                conductivity(point.x, point.y), "Conductivity"
            );
            const double weight = std::abs(map.determinant) * q.weight;
            const auto basis = evaluate_triangle_h1_basis(order, q.r, q.s);
            std::array<std::array<double, 2>, 6> gradients{};
            for (std::size_t i = 0; i < count; ++i) {
                const auto& grad = basis.reference_gradients[i];
                gradients[i] = {
                    (map.j11 * grad[0] - map.j10 * grad[1]) / map.determinant,
                    (-map.j01 * grad[0] + map.j00 * grad[1]) / map.determinant
                };
            }
            for (std::size_t i = 0; i < count; ++i) {
                for (std::size_t j = 0; j < count; ++j) {
                    mass[i][j] += c * basis.values[i] * basis.values[j] * weight;
                    stiffness[i][j] += k * (
                        gradients[i][0] * gradients[j][0] +
                        gradients[i][1] * gradients[j][1]
                    ) * weight;
                }
            }
        }
        for (std::size_t i = 0; i < count; ++i) {
            for (std::size_t j = 0; j < count; ++j) {
                if (!std::isfinite(mass[i][j]) || !std::isfinite(stiffness[i][j])) {
                    throw std::overflow_error("Element matrix accumulation overflowed.");
                }
                mass_entries.push_back({ids[i], ids[j], mass[i][j]});
                stiffness_entries.push_back({ids[i], ids[j], stiffness[i][j]});
            }
        }
    }
    return {
        linalg::CsrMatrix::from_triplets(n, n, mass_entries),
        linalg::CsrMatrix::from_triplets(n, n, stiffness_entries)
    };
}

TransientHeat2DResult solve_transient_heat_2d(
    const TransientHeat2DProblem& problem,
    const TransientHeat2DOptions& options
) {
    if (!std::isfinite(options.time_step) || options.time_step <= 0.0 ||
        options.steps == 0U) {
        throw std::invalid_argument("Transient time step and number of steps must be positive.");
    }
    const double final_time =
        options.time_step * static_cast<double>(options.steps);
    if (!std::isfinite(final_time)) {
        throw std::invalid_argument("Transient final time must be finite.");
    }
    if (!problem.source || !problem.initial_temperature) {
        throw std::invalid_argument("Transient source and initial field are required.");
    }
    const H1FunctionSpace2D space(problem.mesh, options.order);
    const ScalarH1Operators2D operators = assemble_scalar_h1_operators_2d(
        problem.mesh, options.order,
        problem.volumetric_capacity, problem.conductivity
    );
    const BoundaryData boundary = prepare_boundary_data(problem, space);
    const std::size_t n = space.dimension();
    const linalg::CsrMatrix full_system = linalg::CsrMatrix::linear_combination(
        operators.mass, 1.0, operators.stiffness, options.time_step
    );

    std::vector<linalg::Triplet> constrained_entries;
    constrained_entries.reserve(full_system.non_zeros() + n);
    std::vector<linalg::Triplet> boundary_couplings;
    full_system.for_each_nonzero([&](std::size_t row, std::size_t col, double value) {
        if (!boundary.constrained[row] && !boundary.constrained[col]) {
            constrained_entries.push_back({row, col, value});
        } else if (!boundary.constrained[row] && boundary.constrained[col]) {
            boundary_couplings.push_back({row, col, value});
        }
    });
    for (std::size_t i = 0; i < n; ++i) {
        if (boundary.constrained[i]) {
            constrained_entries.push_back({i, i, 1.0});
        }
    }
    const linalg::CsrMatrix system = linalg::CsrMatrix::from_triplets(
        n, n, constrained_entries
    );
    std::vector<double> state = space.interpolate(problem.initial_temperature);
    for (double value : state) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Initial temperature must be finite.");
        }
    }
    const auto initial_boundary = boundary_values(boundary, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (boundary.constrained[i]) {
            state[i] = initial_boundary[i];
        }
    }
    TransientHeat2DResult result;
    result.final_time = final_time;
    result.linear_solvers.reserve(options.steps);
    for (std::size_t step = 1; step <= options.steps; ++step) {
        const double time = static_cast<double>(step) * options.time_step;
        const auto prescribed = boundary_values(boundary, time);
        std::vector<double> rhs = operators.mass.multiply(std::span<const double>(state));
        const auto load = assemble_load(space, problem.source, time);
        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] += options.time_step * load[i];
            if (!std::isfinite(rhs[i])) {
                throw std::overflow_error("Transient right-hand side overflowed.");
            }
        }
        for (const auto& entry : boundary_couplings) {
            rhs[entry.row] -= entry.value * prescribed[entry.column];
            if (!std::isfinite(rhs[entry.row])) {
                throw std::overflow_error("Transient Dirichlet RHS correction overflowed.");
            }
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (boundary.constrained[i]) {
                rhs[i] = prescribed[i];
            }
        }
        linalg::SolverResult step_result = linalg::solve_conjugate_gradient(
            system, std::span<const double>(rhs), options.solver_options
        );
        if (!step_result.converged) {
            throw std::runtime_error(
                "Transient linear solve failed at step " +
                std::to_string(step) + ": " + step_result.message
            );
        }
        state = step_result.solution;
        result.linear_solvers.push_back(std::move(step_result));
    }
    result.temperature = std::move(state);
    return result;
}

} // namespace nova::fem
