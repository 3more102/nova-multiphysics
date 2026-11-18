#include "test_framework.hpp"

#include "nova/fem/transient_heat_2d.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using nova::fem::H1TriangleOrder;
using nova::fem::TransientHeat2DProblem;
using nova::fem::TransientHeat2DOptions;

TransientHeat2DProblem make_problem(std::size_t subdivisions = 3U) {
    auto grid = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, subdivisions, subdivisions
    );
    TransientHeat2DProblem problem;
    problem.mesh = std::move(grid.mesh);
    problem.volumetric_capacity = [](double, double) { return 1.0; };
    problem.conductivity = [](double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    problem.initial_temperature = [](double, double) { return 0.0; };
    return problem;
}

void set_all_dirichlet(
    TransientHeat2DProblem& problem,
    const nova::fem::TransientScalarField2D& temperature
) {
    // Structured boundary edges are inferred from exterior edge incidence.
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> incidence;
    for (const auto& triangle : problem.mesh.elements) {
        const auto add = [&](std::size_t first, std::size_t second) {
            ++incidence[std::minmax(first, second)];
        };
        add(triangle[0], triangle[1]);
        add(triangle[1], triangle[2]);
        add(triangle[2], triangle[0]);
    }
    for (const auto& [edge, count] : incidence) {
        if (count == 1U) {
            problem.dirichlet_edges.push_back({
                {edge.first, edge.second}, temperature
            });
        }
    }
}

double exact_error(
    const TransientHeat2DProblem& problem,
    const nova::fem::TransientHeat2DResult& result,
    H1TriangleOrder order,
    const nova::fem::SpatialScalarField2D& exact
) {
    nova::fem::H1FunctionSpace2D space(problem.mesh, order);
    const auto expected = space.interpolate(exact);
    NOVA_REQUIRE(result.temperature.size() == expected.size());
    double squared = 0.0;
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const double difference = result.temperature[i] - expected[i];
        squared += difference * difference;
    }
    return std::sqrt(squared / static_cast<double>(expected.size()));
}

void check_linear_time_solution(H1TriangleOrder order) {
    auto problem = make_problem();
    // Exact u(x,y,t)=x+y+t, with source rho*c*u_t = 1.
    problem.source = [](double, double, double) { return 1.0; };
    problem.initial_temperature = [](double x, double y) { return x + y; };
    set_all_dirichlet(problem,
        [](double x, double y, double t) { return x + y + t; });

    TransientHeat2DOptions options;
    options.order = order;
    options.steps = 5U;
    options.time_step = 0.2;
    const auto result = nova::fem::solve_transient_heat_2d(problem, options);
    NOVA_REQUIRE_NEAR(result.final_time, 1.0, 1.0e-14);
    NOVA_REQUIRE(result.linear_solvers.size() == options.steps);
    for (const auto& solver : result.linear_solvers) {
        NOVA_REQUIRE(solver.converged);
    }
    NOVA_REQUIRE(
        exact_error(problem, result, order, [](double x, double y) {
            return x + y + 1.0;
        }) < 2.0e-9
    );
}

double quadratic_time_error(std::size_t steps) {
    auto problem = make_problem();
    // Exact P2 spatial polynomial u(x,y,t)=t^2*x*(1-x).
    // f=rho*c*u_t-div(k grad u)=2t*x*(1-x)+2t^2.
    problem.source = [](double x, double, double t) {
        return 2.0 * t * x * (1.0 - x) + 2.0 * t * t;
    };
    set_all_dirichlet(problem, [](double x, double, double t) {
        return t * t * x * (1.0 - x);
    });
    TransientHeat2DOptions options;
    options.order = H1TriangleOrder::P2;
    options.steps = steps;
    options.time_step = 1.0 / static_cast<double>(steps);
    const auto result = nova::fem::solve_transient_heat_2d(problem, options);
    return exact_error(problem, result, options.order,
        [](double x, double) { return x * (1.0 - x); });
}

} // namespace

NOVA_TEST(scalar_h1_p1_mass_matrix_exact_on_reference_triangle) {
    nova::fem::TriangleMesh2D mesh;
    mesh.nodes = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
    mesh.elements = {{0, 1, 2}};
    const auto operators = nova::fem::assemble_scalar_h1_operators_2d(
        mesh, H1TriangleOrder::P1,
        [](double, double) { return 1.0; },
        [](double, double) { return 1.0; }
    );
    NOVA_REQUIRE(operators.mass.rows() == 3U);
    for (std::size_t i = 0; i < 3U; ++i) {
        for (std::size_t j = 0; j < 3U; ++j) {
            const double expected = i == j ? 1.0 / 12.0 : 1.0 / 24.0;
            NOVA_REQUIRE_NEAR(operators.mass.at(i, j), expected, 1.0e-13);
            NOVA_REQUIRE_NEAR(
                operators.mass.at(i, j), operators.mass.at(j, i), 1.0e-14
            );
        }
    }
}

NOVA_TEST(scalar_h1_p2_mass_matrix_preserves_unit_integral) {
    nova::fem::TriangleMesh2D mesh;
    mesh.nodes = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
    mesh.elements = {{0, 1, 2}};
    const auto operators = nova::fem::assemble_scalar_h1_operators_2d(
        mesh, H1TriangleOrder::P2,
        [](double, double) { return 2.0; },
        [](double, double) { return 1.0; }
    );
    NOVA_REQUIRE(operators.mass.rows() == 6U);
    const std::vector<double> ones(6U, 1.0);
    const auto action = operators.mass.multiply(ones);
    double total = 0.0;
    for (std::size_t i = 0; i < 6U; ++i) {
        total += action[i];
        NOVA_REQUIRE(operators.mass.at(i, i) > 0.0);
        for (std::size_t j = 0; j < 6U; ++j) {
            NOVA_REQUIRE_NEAR(
                operators.mass.at(i, j), operators.mass.at(j, i), 1.0e-13
            );
        }
    }
    NOVA_REQUIRE_NEAR(total, 1.0, 1.0e-12); // area 1/2 times capacity 2
    const auto stiffness_action = operators.stiffness.multiply(ones);
    for (const double value : stiffness_action) {
        NOVA_REQUIRE_NEAR(value, 0.0, 1.0e-12);
    }
}

NOVA_TEST(transient_heat_p1_exact_linear_spacetime_solution) {
    check_linear_time_solution(H1TriangleOrder::P1);
}

NOVA_TEST(transient_heat_p2_exact_linear_spacetime_solution) {
    check_linear_time_solution(H1TriangleOrder::P2);
}

NOVA_TEST(transient_heat_p2_first_order_temporal_convergence) {
    const double coarse = quadratic_time_error(4U);
    const double medium = quadratic_time_error(8U);
    const double fine = quadratic_time_error(16U);
    NOVA_REQUIRE(coarse > medium);
    NOVA_REQUIRE(medium > fine);
    NOVA_REQUIRE(coarse / medium > 1.5);
    NOVA_REQUIRE(medium / fine > 1.5);
}

NOVA_TEST(transient_heat_constant_solution_with_zero_flux_boundary) {
    auto problem = make_problem(2U);
    problem.initial_temperature = [](double, double) { return 3.25; };
    for (const H1TriangleOrder order : {H1TriangleOrder::P1, H1TriangleOrder::P2}) {
        TransientHeat2DOptions options;
        options.order = order;
        options.steps = 4U;
        options.time_step = 0.1;
        const auto result = nova::fem::solve_transient_heat_2d(problem, options);
        for (const double value : result.temperature) {
            NOVA_REQUIRE_NEAR(value, 3.25, 1.0e-9);
        }
    }
}

NOVA_TEST(transient_heat_rejects_bad_time_step) {
    auto problem = make_problem();
    TransientHeat2DOptions options;
    options.time_step = 0.0;
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_transient_heat_2d(problem, options));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(transient_heat_rejects_non_positive_capacity) {
    auto problem = make_problem();
    problem.volumetric_capacity = [](double, double) { return -1.0; };
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_transient_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(transient_heat_rejects_interior_dirichlet_edge) {
    auto problem = make_problem(1U);
    problem.dirichlet_edges.push_back({
        {0, 3}, [](double, double, double) { return 0.0; }
    });
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_transient_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(transient_heat_rejects_nonfinite_source) {
    auto problem = make_problem();
    problem.source = [](double, double, double) {
        return std::numeric_limits<double>::infinity();
    };
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_transient_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(transient_heat_rejects_conflicting_corner_temperatures) {
    auto problem = make_problem(1U);
    set_all_dirichlet(problem,
        [](double, double, double) { return 0.0; });
    NOVA_REQUIRE(problem.dirichlet_edges.size() > 2U);
    problem.dirichlet_edges.front().temperature =
        [](double, double, double) { return 1.0; };
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_transient_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}


namespace {

double crank_nicolson_quadratic_time_error(std::size_t steps) {
    auto problem = make_problem(3U);
    // Spatial P2 exact: u=t^2*x*(1-x). Crank-Nicolson is second-order in t.
    problem.source = [](double x, double, double t) {
        return 2.0 * t * x * (1.0 - x) + 2.0 * t * t;
    };
    set_all_dirichlet(problem, [](double x, double, double t) {
        return t * t * x * (1.0 - x);
    });
    TransientHeat2DOptions options;
    options.order = H1TriangleOrder::P2;
    options.scheme = nova::fem::TransientTimeScheme::CrankNicolson;
    options.steps = steps;
    options.time_step = 1.0 / static_cast<double>(steps);
    const auto result = nova::fem::solve_transient_heat_2d(problem, options);
    for (const auto& item : result.diagnostics) {
        NOVA_REQUIRE(std::abs(item.balance_residual) < 2.0e-9);
    }
    return exact_error(problem, result, H1TriangleOrder::P2,
        [](double x, double) { return x * (1.0 - x); });
}

std::vector<nova::fem::BoundaryEdge2D> right_boundary_edges(
    const TransientHeat2DProblem& problem
) {
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> incidence;
    for (const auto& element : problem.mesh.elements) {
        ++incidence[std::minmax(element[0], element[1])];
        ++incidence[std::minmax(element[1], element[2])];
        ++incidence[std::minmax(element[2], element[0])];
    }
    std::vector<nova::fem::BoundaryEdge2D> result;
    for (const auto& [edge, count] : incidence) {
        const auto& a = problem.mesh.nodes[edge.first];
        const auto& b = problem.mesh.nodes[edge.second];
        if (count == 1U && a.x == 1.0 && b.x == 1.0) {
            result.push_back({edge.first, edge.second});
        }
    }
    return result;
}

} // namespace

NOVA_TEST(transient_heat_crank_nicolson_second_order_temporal_convergence) {
    const double coarse = crank_nicolson_quadratic_time_error(4U);
    const double medium = crank_nicolson_quadratic_time_error(8U);
    const double fine = crank_nicolson_quadratic_time_error(16U);
    NOVA_REQUIRE(coarse > medium);
    NOVA_REQUIRE(medium > fine);
    NOVA_REQUIRE(coarse / medium > 3.0);
    NOVA_REQUIRE(medium / fine > 3.0);
}

NOVA_TEST(transient_heat_crank_nicolson_exact_linear_spacetime_solution) {
    auto problem = make_problem();
    problem.source = [](double, double, double) { return 1.0; };
    problem.initial_temperature = [](double x, double y) { return x + y; };
    set_all_dirichlet(problem,
        [](double x, double y, double t) { return x + y + t; });
    for (const H1TriangleOrder order : {H1TriangleOrder::P1, H1TriangleOrder::P2}) {
        TransientHeat2DOptions options;
        options.order = order;
        options.scheme = nova::fem::TransientTimeScheme::CrankNicolson;
        options.steps = 5U;
        options.time_step = 0.2;
        const auto result = nova::fem::solve_transient_heat_2d(problem, options);
        NOVA_REQUIRE(exact_error(problem, result, order,
            [](double x, double y) { return x + y + 1.0; }) < 2.0e-9);
    }
}

NOVA_TEST(transient_heat_source_only_satisfies_global_energy_budget) {
    for (const auto scheme : {
        nova::fem::TransientTimeScheme::BackwardEuler,
        nova::fem::TransientTimeScheme::CrankNicolson
    }) {
        auto problem = make_problem(2U);
        problem.volumetric_capacity = [](double, double) { return 2.0; };
        problem.source = [](double, double, double) { return 2.0; };
        problem.initial_temperature = [](double, double) { return 0.5; };
        TransientHeat2DOptions options;
        options.order = H1TriangleOrder::P2;
        options.scheme = scheme;
        options.steps = 4U;
        options.time_step = 0.1;
        const auto result = nova::fem::solve_transient_heat_2d(problem, options);
        NOVA_REQUIRE(result.diagnostics.size() == options.steps);
        for (const double temperature : result.temperature) {
            NOVA_REQUIRE_NEAR(temperature, 0.9, 1.0e-9);
        }
        for (std::size_t step = 0; step < options.steps; ++step) {
            const auto& d = result.diagnostics[step];
            NOVA_REQUIRE_NEAR(d.heat_content, 1.2 + 0.2 * static_cast<double>(step), 1.0e-9);
            NOVA_REQUIRE_NEAR(d.source_energy, 0.2, 1.0e-10);
            NOVA_REQUIRE_NEAR(d.neumann_energy, 0.0, 1.0e-12);
            NOVA_REQUIRE_NEAR(d.robin_energy, 0.0, 1.0e-12);
            NOVA_REQUIRE_NEAR(d.dirichlet_energy, 0.0, 1.0e-12);
            NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1.0e-9);
        }
    }
}

NOVA_TEST(transient_heat_time_dependent_neumann_flux_energy_budget) {
    for (const auto scheme : {
        nova::fem::TransientTimeScheme::BackwardEuler,
        nova::fem::TransientTimeScheme::CrankNicolson
    }) {
        auto problem = make_problem(2U);
        const auto edges = right_boundary_edges(problem);
        NOVA_REQUIRE(edges.size() == 2U);
        for (const auto& edge : edges) {
            problem.neumann_edges.push_back({
                edge, [](double, double, double t) { return 1.0 + t; }
            });
        }
        TransientHeat2DOptions options;
        options.order = H1TriangleOrder::P2;
        options.scheme = scheme;
        options.steps = 4U;
        options.time_step = 0.25;
        const auto result = nova::fem::solve_transient_heat_2d(problem, options);
        double integrated_flux = 0.0;
        for (std::size_t step = 0; step < options.steps; ++step) {
            const auto& d = result.diagnostics[step];
            const double old_time = 0.25 * static_cast<double>(step);
            const double new_time = old_time + 0.25;
            const double flux_average = scheme == nova::fem::TransientTimeScheme::BackwardEuler
                ? (1.0 + new_time)
                : (1.0 + 0.5 * (old_time + new_time));
            integrated_flux -= options.time_step * flux_average;
            NOVA_REQUIRE_NEAR(d.neumann_energy, -options.time_step * flux_average, 1.0e-10);
            NOVA_REQUIRE_NEAR(d.heat_content, integrated_flux, 1.0e-9);
            NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1.0e-9);
        }
    }
}

NOVA_TEST(transient_heat_time_dependent_robin_flux_energy_budget) {
    auto problem = make_problem(2U);
    for (const auto& edge : right_boundary_edges(problem)) {
        problem.robin_edges.push_back({
            edge,
            [](double, double, double t) { return 2.0 + t; },
            [](double, double, double t) { return 1.0 + t; }
        });
    }
    TransientHeat2DOptions options;
    options.order = H1TriangleOrder::P2;
    options.scheme = nova::fem::TransientTimeScheme::CrankNicolson;
    options.steps = 3U;
    options.time_step = 0.1;
    const auto result = nova::fem::solve_transient_heat_2d(problem, options);
    double accumulated = 0.0;
    for (const auto& d : result.diagnostics) {
        NOVA_REQUIRE(d.robin_energy > 0.0);
        NOVA_REQUIRE_NEAR(d.neumann_energy, 0.0, 1.0e-12);
        accumulated += d.robin_energy;
        NOVA_REQUIRE_NEAR(d.heat_content, accumulated, 2.0e-9);
        NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 2.0e-9);
    }
}

NOVA_TEST(transient_heat_dirichlet_reaction_closes_energy_budget) {
    auto problem = make_problem(2U);
    set_all_dirichlet(problem, [](double, double, double t) { return t; });
    TransientHeat2DOptions options;
    options.scheme = nova::fem::TransientTimeScheme::CrankNicolson;
    options.time_step = 0.1;
    options.steps = 3U;
    const auto result = nova::fem::solve_transient_heat_2d(problem, options);
    for (const auto& d : result.diagnostics) {
        NOVA_REQUIRE(d.dirichlet_energy > 0.0);
        NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 2.0e-9);
    }
}

NOVA_TEST(transient_heat_rejects_dirichlet_neumann_overlap) {
    auto problem = make_problem(1U);
    const auto edge = right_boundary_edges(problem).front();
    problem.dirichlet_edges.push_back({
        edge, [](double, double, double) { return 0.0; }
    });
    problem.neumann_edges.push_back({
        edge, [](double, double, double) { return 1.0; }
    });
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_transient_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(transient_heat_rejects_neumann_robin_overlap) {
    auto problem = make_problem(1U);
    const auto edge = right_boundary_edges(problem).front();
    problem.neumann_edges.push_back({
        edge, [](double, double, double) { return 1.0; }
    });
    problem.robin_edges.push_back({
        edge,
        [](double, double, double) { return 2.0; },
        [](double, double, double) { return 1.0; }
    });
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_transient_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(transient_heat_rejects_negative_robin_coefficient) {
    auto problem = make_problem(1U);
    const auto edge = right_boundary_edges(problem).front();
    problem.robin_edges.push_back({
        edge,
        [](double, double, double) { return -1.0; },
        [](double, double, double) { return 0.0; }
    });
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_transient_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(transient_heat_rejects_unsupported_time_scheme) {
    auto problem = make_problem(1U);
    TransientHeat2DOptions options;
    options.scheme = static_cast<nova::fem::TransientTimeScheme>(999);
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_transient_heat_2d(problem, options));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}
