#include "test_framework.hpp"

#include "nova/fem/nonlinear_heat_2d.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using nova::fem::NonlinearHeat2DProblem;
using nova::fem::NonlinearHeat2DOptions;
using nova::fem::H1TriangleOrder;

NonlinearHeat2DProblem make_nonlinear_problem(std::size_t cells = 3U) {
    auto generated = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, cells, cells
    );
    NonlinearHeat2DProblem problem;
    problem.mesh = std::move(generated.mesh);
    problem.volumetric_capacity = [](double, double) { return 1.0; };
    problem.conductivity = [](double, double, double t) {
        return 1.0 + t * t;
    };
    problem.conductivity_temperature_derivative =
        [](double, double, double t) { return 2.0 * t; };
    problem.source = [](double, double, double) { return 0.0; };
    problem.initial_temperature = [](double, double) { return 0.0; };
    return problem;
}

void apply_all_dirichlet(
    NonlinearHeat2DProblem& problem,
    const nova::fem::TransientScalarField2D& value
) {
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> incidence;
    for (const auto& cell : problem.mesh.elements) {
        ++incidence[std::minmax(cell[0], cell[1])];
        ++incidence[std::minmax(cell[1], cell[2])];
        ++incidence[std::minmax(cell[2], cell[0])];
    }
    for (const auto& [edge, number] : incidence) {
        if (number == 1U) {
            problem.dirichlet_edges.push_back({
                {edge.first, edge.second}, value
            });
        }
    }
}

void check_exact_linear_solution(H1TriangleOrder order) {
    auto problem = make_nonlinear_problem(3U);
    // Exact T=x+y+t, with rho*c=1, k(T)=1+T^2.
    // div(k grad(T))=k'(T) |grad(T)|^2=4T.
    problem.initial_temperature = [](double x, double y) {
        return x + y;
    };
    problem.source = [](double x, double y, double t) {
        return 1.0 - 4.0 * (x + y + t);
    };
    apply_all_dirichlet(problem, [](double x, double y, double t) {
        return x + y + t;
    });

    NonlinearHeat2DOptions opts;
    opts.order = order;
    opts.time_step = 0.2;
    opts.steps = 3U;
    opts.newton_options.residual_tolerance = 1.0e-10;
    const auto result = nova::fem::solve_nonlinear_heat_2d(problem, opts);
    NOVA_REQUIRE_NEAR(result.final_time, 0.6, 1.0e-13);
    NOVA_REQUIRE(result.newton_steps.size() == 3U);
    for (const auto& s : result.newton_steps) {
        NOVA_REQUIRE(s.converged);
        NOVA_REQUIRE(s.residual_norm < 1.0e-9);
    }
    nova::fem::H1FunctionSpace2D space(problem.mesh, order);
    const auto exact = space.interpolate(
        [](double x, double y) { return x + y + 0.6; }
    );
    NOVA_REQUIRE(exact.size() == result.temperature.size());
    for (std::size_t i = 0; i < exact.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.temperature[i], exact[i], 1.0e-8);
    }
}

} // namespace

NOVA_TEST(nonlinear_heat_p1_exact_manufactured_spacetime_solution) {
    check_exact_linear_solution(H1TriangleOrder::P1);
}

NOVA_TEST(nonlinear_heat_p2_exact_manufactured_spacetime_solution) {
    check_exact_linear_solution(H1TriangleOrder::P2);
}

NOVA_TEST(nonlinear_heat_p2_exact_jacobian_matches_finite_difference) {
    auto problem = make_nonlinear_problem(2U);
    apply_all_dirichlet(problem, [](double x, double y, double t) {
        return x + y + t;
    });
    // Include Robin/Neumann boundary data in a separate mixed-edge
    // test below. This test probes free/constrained rows and k'(T).
    const H1TriangleOrder order = H1TriangleOrder::P2;
    nova::fem::H1FunctionSpace2D space(problem.mesh, order);
    const auto old_state = space.interpolate([](double x, double y) {
        return 0.4 + 0.2 * x + 0.3 * y;
    });
    const auto trial = space.interpolate([](double x, double y) {
        return 0.6 + 0.4 * x + 0.25 * y + 0.12 * x * x;
    });
    std::vector<double> direction(trial.size(), 0.0);
    for (std::size_t i = 0; i < direction.size(); ++i) {
        direction[i] = 0.2 + 0.07 * static_cast<double>(i % 5U);
    }
    const auto assembled = nova::fem::assemble_nonlinear_heat_be_system_2d(
        problem, order, old_state, trial, 0.15, 0.15
    );
    const auto jv = assembled.jacobian.multiply(direction);
    const double eps = 1.0e-6;
    std::vector<double> plus = trial;
    std::vector<double> minus = trial;
    for (std::size_t i = 0; i < plus.size(); ++i) {
        plus[i] += eps * direction[i];
        minus[i] -= eps * direction[i];
    }
    const auto r_plus = nova::fem::assemble_nonlinear_heat_be_system_2d(
        problem, order, old_state, plus, 0.15, 0.15
    ).residual;
    const auto r_minus = nova::fem::assemble_nonlinear_heat_be_system_2d(
        problem, order, old_state, minus, 0.15, 0.15
    ).residual;
    for (std::size_t i = 0; i < jv.size(); ++i) {
        const double finite_difference =
            (r_plus[i] - r_minus[i]) / (2.0 * eps);
        NOVA_REQUIRE_NEAR(jv[i], finite_difference, 2.0e-7);
    }
}

NOVA_TEST(nonlinear_heat_p2_accepts_reversed_triangle_orientation) {
    auto problem = make_nonlinear_problem(2U);
    std::swap(problem.mesh.elements[0][1], problem.mesh.elements[0][2]);
    problem.initial_temperature = [](double x, double y) { return x + y; };
    problem.source = [](double x, double y, double t) {
        return 1.0 - 4.0 * (x + y + t);
    };
    apply_all_dirichlet(problem, [](double x, double y, double t) {
        return x + y + t;
    });
    NonlinearHeat2DOptions opts;
    opts.order = H1TriangleOrder::P2;
    opts.time_step = 0.1;
    opts.steps = 2U;
    const auto result = nova::fem::solve_nonlinear_heat_2d(problem, opts);
    nova::fem::H1FunctionSpace2D space(problem.mesh, opts.order);
    const auto expected = space.interpolate([](double x, double y) {
        return x + y + 0.2;
    });
    for (std::size_t i = 0; i < expected.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.temperature[i], expected[i], 1.0e-8);
    }
}

NOVA_TEST(nonlinear_heat_homogeneous_neumann_constant_heating) {
    auto problem = make_nonlinear_problem(2U);
    problem.initial_temperature = [](double, double) { return 0.5; };
    problem.source = [](double, double, double) { return 2.0; };
    NonlinearHeat2DOptions opts;
    opts.order = H1TriangleOrder::P2;
    opts.steps = 3U;
    opts.time_step = 0.1;
    const auto result = nova::fem::solve_nonlinear_heat_2d(problem, opts);
    for (const double value : result.temperature) {
        NOVA_REQUIRE_NEAR(value, 1.1, 1.0e-8);
    }
}

NOVA_TEST(nonlinear_heat_p1_robin_boundary_with_constant_solution) {
    auto problem = make_nonlinear_problem(2U);
    problem.initial_temperature = [](double, double) { return 2.0; };
    // Matching ambient makes Robin flux identically zero.
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> incidence;
    for (const auto& tri : problem.mesh.elements) {
        ++incidence[std::minmax(tri[0], tri[1])];
        ++incidence[std::minmax(tri[1], tri[2])];
        ++incidence[std::minmax(tri[2], tri[0])];
    }
    for (const auto& [edge, count] : incidence) {
        if (count == 1U) {
            problem.robin_edges.push_back({
                {edge.first, edge.second},
                [](double, double, double) { return 5.0; },
                [](double, double, double) { return 2.0; }
            });
        }
    }
    NonlinearHeat2DOptions opts;
    opts.order = H1TriangleOrder::P1;
    opts.time_step = 0.2;
    const auto result = nova::fem::solve_nonlinear_heat_2d(problem, opts);
    for (const double value : result.temperature) {
        NOVA_REQUIRE_NEAR(value, 2.0, 1.0e-9);
    }
}

NOVA_TEST(nonlinear_heat_rejects_negative_conductivity) {
    auto problem = make_nonlinear_problem(1U);
    problem.conductivity = [](double, double, double) { return -0.5; };
    const auto state = nova::fem::H1FunctionSpace2D(
        problem.mesh, H1TriangleOrder::P1
    ).interpolate([](double, double) { return 0.0; });
    bool caught = false;
    try {
        static_cast<void>(nova::fem::assemble_nonlinear_heat_be_system_2d(
            problem, H1TriangleOrder::P1, state, state, 0.1, 0.1
        ));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(nonlinear_heat_rejects_nonfinite_conductivity_derivative) {
    auto problem = make_nonlinear_problem(1U);
    problem.conductivity_temperature_derivative =
        [](double, double, double) {
            return std::numeric_limits<double>::quiet_NaN();
        };
    const auto state = nova::fem::H1FunctionSpace2D(
        problem.mesh, H1TriangleOrder::P1
    ).interpolate([](double, double) { return 0.0; });
    bool caught = false;
    try {
        static_cast<void>(nova::fem::assemble_nonlinear_heat_be_system_2d(
            problem, H1TriangleOrder::P1, state, state, 0.1, 0.1
        ));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(nonlinear_heat_rejects_invalid_dof_dimensions) {
    auto problem = make_nonlinear_problem(1U);
    const std::vector<double> wrong(1U, 0.0);
    bool caught = false;
    try {
        static_cast<void>(nova::fem::assemble_nonlinear_heat_be_system_2d(
            problem, H1TriangleOrder::P2, wrong, wrong, 0.1, 0.1
        ));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(nonlinear_heat_fails_loudly_on_newton_iteration_limit) {
    auto problem = make_nonlinear_problem(3U);
    problem.initial_temperature = [](double x, double y) { return x + y; };
    problem.source = [](double x, double y, double t) {
        return 1.0 - 4.0 * (x + y + t);
    };
    apply_all_dirichlet(problem, [](double x, double y, double t) {
        return x + y + t;
    });
    NonlinearHeat2DOptions opts;
    opts.order = H1TriangleOrder::P2;
    opts.time_step = 0.4;
    opts.newton_options.max_iterations = 1U;
    opts.newton_options.residual_tolerance = 1.0e-15;
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_nonlinear_heat_2d(problem, opts));
    } catch (const std::runtime_error&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}


namespace {

void check_nonlinear_mixed_boundary(
    H1TriangleOrder order, bool use_robin
) {
    auto problem = make_nonlinear_problem(3U);
    problem.initial_temperature = [](double x, double y) { return x + y; };
    problem.source = [](double x, double y, double t) {
        return 1.0 - 4.0 * (x + y + t);
    };
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> incidence;
    for (const auto& cell : problem.mesh.elements) {
        ++incidence[std::minmax(cell[0], cell[1])];
        ++incidence[std::minmax(cell[1], cell[2])];
        ++incidence[std::minmax(cell[2], cell[0])];
    }
    for (const auto& [edge, count] : incidence) {
        if (count != 1U) {
            continue;
        }
        const auto& a = problem.mesh.nodes[edge.first];
        const auto& b = problem.mesh.nodes[edge.second];
        const nova::fem::BoundaryEdge2D boundary{{
            edge.first, edge.second
        }};
        if (a.x == 1.0 && b.x == 1.0) {
            if (use_robin) {
                // On x=1, n=(1,0). Exact outward q=-k(T).
                // Choose h=2 and ambient=T+k(T)/2.
                problem.robin_edges.push_back({
                    boundary,
                    [](double, double, double) { return 2.0; },
                    [](double x, double y, double time) {
                        const double t = x + y + time;
                        return t + 0.5 * (1.0 + t * t);
                    }
                });
            } else {
                problem.neumann_edges.push_back({
                    boundary,
                    [](double x, double y, double time) {
                        const double t = x + y + time;
                        return -(1.0 + t * t);
                    }
                });
            }
        } else {
            problem.dirichlet_edges.push_back({
                boundary,
                [](double x, double y, double time) {
                    return x + y + time;
                }
            });
        }
    }

    NonlinearHeat2DOptions options;
    options.order = order;
    options.time_step = 0.1;
    options.steps = 2U;
    const auto result = nova::fem::solve_nonlinear_heat_2d(problem, options);
    nova::fem::H1FunctionSpace2D space(problem.mesh, order);
    const auto exact = space.interpolate([](double x, double y) {
        return x + y + 0.2;
    });
    NOVA_REQUIRE(exact.size() == result.temperature.size());
    for (std::size_t i = 0; i < exact.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.temperature[i], exact[i], 3.0e-8);
    }
}

} // namespace

NOVA_TEST(nonlinear_heat_p1_p2_exact_manufactured_neumann_flux) {
    check_nonlinear_mixed_boundary(H1TriangleOrder::P1, false);
    check_nonlinear_mixed_boundary(H1TriangleOrder::P2, false);
}

NOVA_TEST(nonlinear_heat_p1_p2_exact_manufactured_robin_convection) {
    check_nonlinear_mixed_boundary(H1TriangleOrder::P1, true);
    check_nonlinear_mixed_boundary(H1TriangleOrder::P2, true);
}

NOVA_TEST(nonlinear_heat_rejects_overlapping_boundary_conditions) {
    auto problem = make_nonlinear_problem(1U);
    apply_all_dirichlet(problem, [](double, double, double) {
        return 0.0;
    });
    const auto edge = problem.dirichlet_edges.front().edge;
    problem.neumann_edges.push_back({
        edge, [](double, double, double) { return 0.0; }
    });
    bool rejected = false;
    try {
        static_cast<void>(nova::fem::solve_nonlinear_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}
