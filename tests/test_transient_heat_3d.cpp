#include "test_framework.hpp"
#include "nova/fem/simplex_quadrature.hpp"
#include "nova/fem/transient_heat_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
using namespace nova::fem;

double one(double, double, double) { return 1.0; }
double zero(double, double, double, double) { return 0.0; }

TransientHeat3DProblem reference_tet() {
    TransientHeat3DProblem p;
    p.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    p.volumetric_capacity = one;
    p.conductivity = one;
    p.source = zero;
    p.initial_temperature = [](double x, double, double) { return x; };
    return p;
}

template<class Exception, class Function>
bool throws_as(Function&& action) {
    try {
        action();
    } catch (const Exception&) {
        return true;
    }
    return false;
}

double exact(double x, double y, double z, double time) {
    return 2.0 + x + 2.0 * y - z + time;
}

bool on_plane(const TetrahedralMesh3D& mesh, BoundaryFace3D face,
              double Point3D::*coordinate, double value) {
    return std::all_of(face.begin(), face.end(), [&](std::size_t i) {
        return mesh.nodes[i].*coordinate == value;
    });
}

TransientHeat3DProblem mixed_cube() {
    auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 2, 2, 2);
    TransientHeat3DProblem p;
    p.mesh = box.mesh;
    for (auto& point : p.mesh.nodes) {
        if (point.x == 0.5 && point.y == 0.5 && point.z == 0.5) {
            point = {0.57, 0.38, 0.59};
        }
    }
    for (std::size_t i = 0; i < p.mesh.elements.size(); i += 3) {
        std::swap(p.mesh.elements[i][1], p.mesh.elements[i][2]);
    }
    p.volumetric_capacity = [](double x, double, double) { return 1.0 + x; };
    p.conductivity = [](double x, double, double) { return 2.0 + x; };
    p.source = [](double x, double, double, double) { return x; };
    // The start-time Dirichlet projection corrects the deliberate offset at x=0.
    p.initial_temperature = [](double x, double y, double z) {
        return exact(x, y, z, 0.0) + (x == 0.0 ? 50.0 : 0.0);
    };
    for (auto face : box.boundary_faces) {
        std::swap(face[0], face[2]);
        if (on_plane(p.mesh, face, &Point3D::x, 0.0)) {
            p.dirichlet_faces.push_back({face, exact});
        } else if (on_plane(p.mesh, face, &Point3D::x, 1.0)) {
            double y_sum = 0.0;
            for (auto id : face) {
                y_sum += p.mesh.nodes[id].y;
            }
            if (y_sum < 1.5) {
                p.neumann_faces.push_back({face, [](double, double, double, double) { return -3.0; }});
            } else {
                p.robin_faces.push_back({face,
                    [](double, double, double, double t) { return 2.0 + t; },
                    [](double x, double y, double z, double t) {
                        return exact(x, y, z, t) + 3.0 / (2.0 + t);
                    }});
            }
        } else {
            double slope = 1.0; // z=1 has q_n = +k.
            if (on_plane(p.mesh, face, &Point3D::y, 0.0)) slope = 2.0;
            if (on_plane(p.mesh, face, &Point3D::y, 1.0)) slope = -2.0;
            if (on_plane(p.mesh, face, &Point3D::z, 0.0)) slope = -1.0;
            p.neumann_faces.push_back({face, [slope](double x, double, double, double) {
                return slope * (2.0 + x);
            }});
        }
    }
    return p;
}
} // namespace

NOVA_TEST(tetrahedron_and_face_rules_integrate_all_barycentric_quadratics) {
    // On a simplex of dimension d, the integral of L_i L_j is
    // measure * (1+delta_ij)/((d+1)(d+2)).
    const auto check = [](const auto& rule, std::size_t n, double measure) {
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                double value = 0.0;
                for (const auto& q : rule) {
                    NOVA_REQUIRE(q.weight > 0.0);
                    value += q.weight * q.shape[i] * q.shape[j];
                }
                const double divisor = static_cast<double>(n * (n + 1));
                NOVA_REQUIRE_NEAR(value, measure * (i == j ? 2.0 : 1.0) / divisor, 2.0e-16);
            }
        }
    };
    check(tetrahedron_degree_two_quadrature, 4, 1.0 / 6.0);
    check(triangle_face_degree_two_quadrature, 3, 0.5);
}

NOVA_TEST(heat_3d_mass_and_stiffness_match_exact_tetrahedron_matrices) {
    auto p = reference_tet();
    p.mesh.nodes = {{0,0,0}, {2,0,0}, {0,3,0}, {0,0,4}};
    const std::array<std::array<double, 3>, 4> grad{{
        {{-0.5,-1.0/3.0,-0.25}}, {{0.5,0,0}}, {{0,1.0/3.0,0}}, {{0,0,0.25}}
    }};
    for (int orientation = 0; orientation < 2; ++orientation) {
        const auto a = assemble_scalar_h1_operators_3d(p.mesh,
            [](double, double, double) { return 3.0; },
            [](double x, double y, double z) { return 1.0+x+2.0*y+3.0*z; });
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = 0; j < 4; ++j) {
                double dot = 0.0;
                for (std::size_t d = 0; d < 3; ++d) dot += grad[i][d] * grad[j][d];
                // Volume 4, mean k = 6, constant capacity 3.
                NOVA_REQUIRE_NEAR(a.mass.at(i,j), 12.0 * (i == j ? 2.0 : 1.0) / 20.0, 2e-14);
                NOVA_REQUIRE_NEAR(a.stiffness.at(i,j), 24.0 * dot, 2e-13);
            }
        }
        std::swap(p.mesh.elements[0][1], p.mesh.elements[0][2]);
    }
}

NOVA_TEST(heat_3d_resolves_mixed_faces_and_independent_physical_heat_budget) {
    const auto p = mixed_cube();
    for (auto scheme : {TransientTimeScheme::BackwardEuler, TransientTimeScheme::CrankNicolson}) {
        TransientHeat3DOptions o;
        o.scheme = scheme;
        o.time_step = 0.05;
        o.steps = 4;
        o.solver_options.relative_tolerance = 1.0e-13;
        const auto r = solve_transient_heat_3d(p, o);
        NOVA_REQUIRE(r.linear_solvers.size() == o.steps && r.diagnostics.size() == o.steps);
        for (std::size_t i = 0; i < r.temperature.size(); ++i) {
            const auto& n = p.mesh.nodes[i];
            NOVA_REQUIRE_NEAR(r.temperature[i], exact(n.x,n.y,n.z,r.final_time), 2.0e-10);
        }
        for (const auto& d : r.diagnostics) {
            // Physical integrals from the analytical gradient and the actual face areas.
            NOVA_REQUIRE_NEAR(d.heat_content, 55.0/12.0 + 1.5*d.time, 1.0e-10);
            NOVA_REQUIRE_NEAR(d.source_energy, 0.5 * o.time_step, 1.0e-12);
            NOVA_REQUIRE_NEAR(d.neumann_energy, 1.5 * o.time_step, 1.0e-12);
            NOVA_REQUIRE_NEAR(d.robin_energy, 1.5 * o.time_step, 1.0e-10);
            NOVA_REQUIRE_NEAR(d.dirichlet_energy, -2.0 * o.time_step, 1.0e-10);
            NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1.0e-10);
        }
    }
}

NOVA_TEST(heat_3d_time_orders_match_a_known_semidiscrete_decay_mode) {
    auto p = reference_tet();
    p.dirichlet_faces.push_back({{0,2,3}, zero});
    for (auto scheme : {TransientTimeScheme::BackwardEuler, TransientTimeScheme::CrankNicolson}) {
        double previous_error = 0.0;
        for (auto steps : {10U, 20U, 40U}) {
            TransientHeat3DOptions o;
            o.scheme = scheme;
            o.steps = steps;
            o.time_step = 0.2 / static_cast<double>(steps);
            const auto r = solve_transient_heat_3d(p, o);
            const double z = 10.0 * o.time_step; // K11/M11 = (1/6)/(1/60).
            const double factor = scheme == TransientTimeScheme::BackwardEuler ?
                1.0 / (1.0 + z) : (1.0 - 0.5*z) / (1.0 + 0.5*z);
            NOVA_REQUIRE_NEAR(r.temperature[1], std::pow(factor, steps), 1.0e-13);
            const double error = std::abs(r.temperature[1] - std::exp(-2.0));
            if (previous_error > 0) {
                const double expected = scheme == TransientTimeScheme::BackwardEuler ? 1.0 : 2.0;
                NOVA_REQUIRE_NEAR(std::log2(previous_error / error), expected, 0.08);
            }
            previous_error = error;
        }
    }
}

NOVA_TEST(heat_3d_insulated_components_conserve_their_own_heat_and_assemble_materials_once) {
    auto p = reference_tet();
    for (const auto& point : reference_tet().mesh.nodes) {
        p.mesh.nodes.push_back({point.x + 2.0, point.y, point.z});
    }
    p.mesh.elements.push_back({4,5,6,7});
    p.initial_temperature = [](double x, double, double) { return x < 2.0 ? 3.0 : 5.0; };
    int material_samples = 0;
    p.volumetric_capacity = [&](double, double, double) { ++material_samples; return 1.0; };
    TransientHeat3DOptions o;
    o.steps = 5;
    const auto r = solve_transient_heat_3d(p, o);
    NOVA_REQUIRE(material_samples == 8);
    for (std::size_t i = 0; i < 8; ++i) {
        NOVA_REQUIRE_NEAR(r.temperature[i], i < 4 ? 3.0 : 5.0, 1.0e-12);
    }
    for (const auto& d : r.diagnostics) {
        NOVA_REQUIRE_NEAR(d.heat_content, 8.0/6.0, 1.0e-12);
        NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1.0e-12);
    }
}

NOVA_TEST(heat_3d_restart_reproduces_the_continuous_run_with_time_dependent_robin_data) {
    for (auto scheme : {TransientTimeScheme::BackwardEuler, TransientTimeScheme::CrankNicolson}) {
        auto p = mixed_cube();
        TransientHeat3DOptions o;
        o.scheme = scheme;
        o.steps = 4;
        o.time_step = 0.125;
        const auto full = solve_transient_heat_3d(p, o);
        o.steps = 2;
        const auto first = solve_transient_heat_3d(p, o);
        o.start_time = first.final_time;
        p.initial_dof_values = first.temperature;
        p.initial_temperature = {}; // restart must take precedence
        const auto second = solve_transient_heat_3d(p, o);
        NOVA_REQUIRE(second.final_time == full.final_time);
        for (std::size_t i = 0; i < full.temperature.size(); ++i) {
            NOVA_REQUIRE_NEAR(second.temperature[i], full.temperature[i], 1.0e-12);
        }
    }
}

NOVA_TEST(heat_3d_backward_euler_skips_zero_weight_start_time_volume_and_natural_data) {
    auto p = reference_tet();
    const auto positive_time_zero = [](double, double, double, double time) {
        if (time == 0) throw std::domain_error("No initial endpoint data.");
        return 0.0;
    };
    p.source = positive_time_zero;
    p.neumann_faces.push_back({{0,1,2}, positive_time_zero});
    p.robin_faces.push_back({{1,2,3}, positive_time_zero, positive_time_zero});
    TransientHeat3DOptions o;
    NOVA_REQUIRE(solve_transient_heat_3d(p, o).linear_solvers[0].converged);
    o.scheme = TransientTimeScheme::CrankNicolson;
    NOVA_REQUIRE(throws_as<std::domain_error>([&] { (void)solve_transient_heat_3d(p, o); }));
}

NOVA_TEST(heat_3d_rejects_invalid_coefficients_faces_states_and_time) {
    const auto rejects = [](const TransientHeat3DProblem& p) {
        return throws_as<std::invalid_argument>([&] { (void)solve_transient_heat_3d(p); });
    };
    auto p = reference_tet();
    p.volumetric_capacity = [](double, double, double) { return 0.0; };
    NOVA_REQUIRE(rejects(p));
    p = reference_tet();
    p.conductivity = [](double, double, double) { return std::numeric_limits<double>::quiet_NaN(); };
    NOVA_REQUIRE(rejects(p));
    p = reference_tet();
    p.initial_dof_values = {0,1};
    NOVA_REQUIRE(rejects(p));
    p.initial_dof_values = {0,1,0,std::numeric_limits<double>::infinity()};
    NOVA_REQUIRE(rejects(p));
    p = reference_tet();
    p.dirichlet_faces.push_back({{0,1,2}, zero});
    p.neumann_faces.push_back({{2,0,1}, zero});
    NOVA_REQUIRE(rejects(p));
    p.neumann_faces.clear();
    p.dirichlet_faces.push_back({{0,1,3}, [](double, double, double, double t) { return t; }});
    NOVA_REQUIRE(rejects(p));
    p = mixed_cube();
    p.neumann_faces.push_back({{0,1,13}, zero}); // Includes the center vertex.
    NOVA_REQUIRE(rejects(p));
    p = reference_tet();
    p.robin_faces.push_back({{0,1,2}, [](double, double, double, double) { return -1.0; }, zero});
    NOVA_REQUIRE(rejects(p));
    p = reference_tet();
    TransientHeat3DOptions o;
    o.start_time = 1.0e16;
    o.time_step = 1.0;
    o.steps = 4;
    NOVA_REQUIRE(throws_as<std::invalid_argument>([&] { (void)solve_transient_heat_3d(p, o); }));
}

NOVA_TEST(heat_3d_reports_a_failed_linear_step) {
    auto p = mixed_cube();
    TransientHeat3DOptions o;
    o.solver_options.max_iterations = 1;
    NOVA_REQUIRE(throws_as<std::runtime_error>([&] { (void)solve_transient_heat_3d(p, o); }));
}
