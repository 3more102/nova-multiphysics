#include "nova/fem/transient_heat_2d.hpp"

#include "nova/fem/triangle_p1.hpp"
#include "heat_time_integration.hpp"

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

void validate_natural_boundaries(
    const TransientHeat2DProblem& problem
) {
    const std::size_t count = problem.mesh.nodes.size();
    std::map<Edge, std::size_t> incidence;
    for (const Triangle2D& cell : problem.mesh.elements) {
        ++incidence[edge_key(cell[0], cell[1])];
        ++incidence[edge_key(cell[1], cell[2])];
        ++incidence[edge_key(cell[2], cell[0])];
    }
    std::set<Edge> used;
    const auto validate = [&](BoundaryEdge2D edge) {
        if (edge[0] >= count || edge[1] >= count) {
            throw std::out_of_range("Thermal boundary edge references a missing vertex.");
        }
        if (edge[0] == edge[1]) {
            throw std::invalid_argument("Thermal boundary edge endpoints must differ.");
        }
        const Edge key = edge_key(edge[0], edge[1]);
        const auto it = incidence.find(key);
        if (it == incidence.end() || it->second != 1U) {
            throw std::invalid_argument("Thermal boundary data must reference an exterior edge.");
        }
        if (!used.insert(key).second) {
            throw std::invalid_argument("Thermal edge has duplicate/conflicting boundary conditions.");
        }
    };
    for (const auto& bc : problem.dirichlet_edges) {
        validate(bc.edge);
    }
    for (const auto& bc : problem.neumann_edges) {
        validate(bc.edge);
        if (!bc.outward_flux) {
            throw std::invalid_argument("Neumann outward flux function is required.");
        }
    }
    for (const auto& bc : problem.robin_edges) {
        validate(bc.edge);
        if (!bc.transfer_coefficient || !bc.ambient_temperature) {
            throw std::invalid_argument("Robin transfer and ambient functions are required.");
        }
    }
}

[[nodiscard]] detail::HeatLoads assemble_heat_loads(
    const TransientHeat2DProblem& problem,
    const H1FunctionSpace2D& space,
    const double time
) {
    const std::size_t n = space.dimension();
    std::vector<linalg::Triplet> robin_entries;
    std::vector<double> neumann_rhs(n, 0.0);
    std::vector<double> robin_ambient_rhs(n, 0.0);
    constexpr std::array<double, 3> points{{
        0.5 * (1.0 - 0.77459666924148337704),
        0.5,
        0.5 * (1.0 + 0.77459666924148337704)
    }};
    constexpr std::array<double, 3> weights{{5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0}};
    const std::size_t edge_count =
        space.order() == H1TriangleOrder::P2 ? 3U : 2U;
    const auto edge_data = [&](BoundaryEdge2D edge) {
        const Point2D& a = problem.mesh.nodes[edge[0]];
        const Point2D& b = problem.mesh.nodes[edge[1]];
        const double length = std::hypot(b.x - a.x, b.y - a.y);
        if (!std::isfinite(length) || length <= 0.0) {
            throw std::invalid_argument("Thermal boundary edge has invalid length.");
        }
        const std::array<std::size_t, 3> ids{{
            edge[0], edge[1],
            space.order() == H1TriangleOrder::P2
                ? space.dofs().edge_dof(edge[0], edge[1]).value() : 0U
        }};
        return std::pair{ids, length};
    };
    const auto edge_basis = [&](double t) {
        if (space.order() == H1TriangleOrder::P1) {
            return std::array<double, 3>{{1.0 - t, t, 0.0}};
        }
        return std::array<double, 3>{{
            (1.0 - t) * (1.0 - 2.0 * t),
            t * (2.0 * t - 1.0),
            4.0 * t * (1.0 - t)
        }};
    };
    if (problem.robin_edges.size() > robin_entries.max_size() / (3U * 9U)) {
        throw std::length_error("Thermal Robin matrix is too large.");
    }
    robin_entries.reserve(problem.robin_edges.size() * 3U * 9U);

    for (const auto& bc : problem.neumann_edges) {
        const auto [ids, length] = edge_data(bc.edge);
        const Point2D& a = problem.mesh.nodes[bc.edge[0]];
        const Point2D& b = problem.mesh.nodes[bc.edge[1]];
        for (std::size_t q = 0; q < points.size(); ++q) {
            const double t = points[q];
            const double x = (1.0 - t) * a.x + t * b.x;
            const double y = (1.0 - t) * a.y + t * b.y;
            const double flux = bc.outward_flux(x, y, time);
            if (!std::isfinite(flux)) {
                throw std::invalid_argument("Neumann outward flux must be finite.");
            }
            const auto phi = edge_basis(t);
            for (std::size_t i = 0; i < edge_count; ++i) {
                // Positive outward heat flux removes heat from the domain.
                neumann_rhs[ids[i]] -= flux * phi[i] * length * weights[q];
                if (!std::isfinite(neumann_rhs[ids[i]])) {
                    throw std::overflow_error("Neumann transient load overflowed.");
                }
            }
        }
    }

    for (const auto& bc : problem.robin_edges) {
        const auto [ids, length] = edge_data(bc.edge);
        const Point2D& a = problem.mesh.nodes[bc.edge[0]];
        const Point2D& b = problem.mesh.nodes[bc.edge[1]];
        for (std::size_t q = 0; q < points.size(); ++q) {
            const double t = points[q];
            const double x = (1.0 - t) * a.x + t * b.x;
            const double y = (1.0 - t) * a.y + t * b.y;
            const double h = bc.transfer_coefficient(x, y, time);
            const double ambient = bc.ambient_temperature(x, y, time);
            if (!std::isfinite(h) || h < 0.0) {
                throw std::invalid_argument("Robin transfer coefficient must be finite and nonnegative.");
            }
            if (!std::isfinite(ambient)) {
                throw std::invalid_argument("Robin ambient temperature must be finite.");
            }
            const auto phi = edge_basis(t);
            const double weight = length * weights[q];
            for (std::size_t i = 0; i < edge_count; ++i) {
                robin_ambient_rhs[ids[i]] += h * ambient * phi[i] * weight;
                if (!std::isfinite(robin_ambient_rhs[ids[i]])) {
                    throw std::overflow_error("Robin transient load overflowed.");
                }
                for (std::size_t j = 0; j < edge_count; ++j) {
                    robin_entries.push_back({
                        ids[i], ids[j], h * phi[i] * phi[j] * weight
                    });
                }
            }
        }
    }

    return {
        linalg::CsrMatrix::from_triplets(n, n, robin_entries),
        std::move(neumann_rhs),
        std::move(robin_ambient_rhs),
        assemble_load(space, problem.source, time)
    };
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
    const TransientHeatTimeOptions time_options{
        options.scheme, options.time_step, options.start_time, options.steps, options.solver_options};
    detail::validate_heat_time_grid(time_options);
    if (!problem.source || (problem.initial_dof_values.empty() && !problem.initial_temperature)) {
        throw std::invalid_argument("Transient source and initial field/DOFs are required.");
    }
    const H1FunctionSpace2D space(problem.mesh, options.order);
    const auto boundary = prepare_boundary_data(problem, space);
    validate_natural_boundaries(problem);
    const auto operators = assemble_scalar_h1_operators_2d(
        problem.mesh, options.order, problem.volumetric_capacity, problem.conductivity);
    auto state = problem.initial_dof_values.empty() ? space.interpolate(problem.initial_temperature)
                                                    : problem.initial_dof_values;
    return detail::integrate_heat(
        operators, time_options, std::move(state), boundary.constrained,
        [&](double time) { return boundary_values(boundary, time); },
        [&](double time) { return assemble_heat_loads(problem, space, time); });
}

} // namespace nova::fem
