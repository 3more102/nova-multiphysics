#include "test_framework.hpp"
#include "nova/fem/transient_heat_3d.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace nova::fem;

TransientHeat3DProblem unit_heat() {
    TransientHeat3DProblem problem;
    auto& p = problem.diffusion;
    p.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    p.cell_diffusivity = {1.0};
    p.source = [](double, double, double) { return 0.0; };
    problem.cell_volumetric_capacity = {1.0};
    problem.initial_temperature = {0.0, 0.0, 0.0, 0.0};
    return problem;
}

template<class F>
bool rejects(F&& fn) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}
} // namespace

NOVA_TEST(transient_heat_3d_insulated_uniform_source_ramps_exactly) {
    auto p = unit_heat();
    p.cell_volumetric_capacity = {4.0};
    p.diffusion.source = [](double, double, double) { return 8.0; };
    p.initial_temperature = {300,300,300,300};
    TransientHeat3DOptions options;
    options.time_step = 0.125;
    options.steps = 10;
    options.solver_options.relative_tolerance = 1e-13;
    options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                               TransientHeat3DScheme::CrankNicolson}) {
        options.scheme = scheme;
        const auto result = solve_transient_heat_3d(p, options);
        NOVA_REQUIRE(result.linear_solvers.size() == options.steps);
        NOVA_REQUIRE(result.diagnostics.size() == options.steps);
        NOVA_REQUIRE_NEAR(result.final_time, 1.25, 1e-14);
        for (double temp : result.temperature) {
            NOVA_REQUIRE_NEAR(temp, 302.5, 2e-9);
        }
        for (const auto& diag : result.diagnostics) {
            NOVA_REQUIRE_NEAR(diag.source_energy, 1.0 / 6.0, 1e-13);
            NOVA_REQUIRE_NEAR(diag.neumann_energy, 0.0, 1e-14);
            NOVA_REQUIRE_NEAR(diag.dirichlet_energy, 0.0, 1e-14);
            NOVA_REQUIRE_NEAR(diag.balance_residual, 0.0, 1e-10);
        }
    }
}

NOVA_TEST(transient_heat_3d_dirichlet_decay_matches_closed_form_time_step) {
    auto p = unit_heat();
    p.initial_temperature = {0,0,0,10};
    p.diffusion.dirichlet_constraints = {{0,0}, {1,0}, {2,0}};
    TransientHeat3DOptions options;
    options.time_step = 0.1;
    options.solver_options.relative_tolerance = 1e-13;
    options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                               TransientHeat3DScheme::CrankNicolson}) {
        options.scheme = scheme;
        const auto result = solve_transient_heat_3d(p, options);
        // For a unit tetra: M_33=V/10=1/60, K_33=V=1/6.
        // BE: 10*(M_33)/(M_33+dt*K_33) = 5.
        // CN: 10*(M_33-dt*K_33/2)/(M_33+dt*K_33/2) = 10/3.
        NOVA_REQUIRE_NEAR(result.temperature[3],
            scheme == TransientHeat3DScheme::BackwardEuler ? 5.0 : 10.0 / 3.0,
            2e-12);
        for (std::size_t i = 0; i < 3; ++i) {
            NOVA_REQUIRE_NEAR(result.temperature[i], 0.0, 1e-15);
        }
        const auto& d = result.diagnostics[0];
        NOVA_REQUIRE(d.dirichlet_energy < 0.0);
        NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 2e-12);
    }
}

NOVA_TEST(transient_heat_3d_neumann_and_robin_satisfy_discrete_energy_budget) {
    auto p = unit_heat();
    p.initial_temperature = {300,300,300,300};
    p.diffusion.neumann_faces.push_back({
        {0,1,2}, [](double, double, double) { return 2.0; }
    });
    TransientHeat3DOptions options;
    options.time_step = 0.1;
    options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    const auto neumann = solve_transient_heat_3d(p, options);
    NOVA_REQUIRE_NEAR(neumann.diagnostics[0].neumann_energy, -0.1, 1e-13);
    NOVA_REQUIRE_NEAR(neumann.diagnostics[0].heat_content, 50.0 - 0.1, 1e-10);
    NOVA_REQUIRE_NEAR(neumann.diagnostics[0].balance_residual, 0.0, 1e-10);

    p.diffusion.neumann_faces.clear();
    p.diffusion.robin_faces.push_back({
        {0,1,2},
        [](double, double, double) { return 3.0; },
        [](double, double, double) { return 290.0; }
    });
    options.scheme = TransientHeat3DScheme::CrankNicolson;
    const auto robin = solve_transient_heat_3d(p, options);
    NOVA_REQUIRE(robin.diagnostics[0].robin_energy < 0.0);
    NOVA_REQUIRE(robin.diagnostics[0].heat_content < 50.0);
    NOVA_REQUIRE_NEAR(robin.diagnostics[0].balance_residual, 0.0, 1e-10);
}

NOVA_TEST(transient_heat_3d_multicell_anisotropic_source_conserves_heat) {
    auto p = unit_heat();
    const auto box = make_box_mesh_3d(0,1,0,1,0,1,1,1,1);
    p.diffusion.mesh = box.mesh;
    p.diffusion.cell_diffusivity.clear();
    p.diffusion.diffusivity_tensor = [](double, double, double) {
        return SymmetricTensor3D{4,3,2,0.3,0.1,-0.2};
    };
    p.diffusion.source = [](double, double, double) { return 1.0; };
    p.diffusion.cell_source.assign(box.mesh.elements.size(), 2.0);
    p.cell_volumetric_capacity.clear();
    p.volumetric_capacity = [](double x, double, double) { return 2.0 + x; };
    p.initial_temperature.assign(box.mesh.nodes.size(), 0.0);
    TransientHeat3DOptions options;
    options.time_step = 0.2;
    options.steps = 3;
    options.solver_options.relative_tolerance = 1e-13;
    options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    const auto result = solve_transient_heat_3d(p, options);
    NOVA_REQUIRE(result.temperature.size() == box.mesh.nodes.size());
    for (const auto& d : result.diagnostics) {
        NOVA_REQUIRE_NEAR(d.source_energy, 0.6, 1e-12);
        NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1e-9);
    }
    NOVA_REQUIRE_NEAR(result.diagnostics.back().heat_content, 1.8, 1e-9);
}

NOVA_TEST(transient_heat_3d_rejects_invalid_heat_capacity_bcs_and_time) {
    auto p = unit_heat();
    NOVA_REQUIRE(rejects([&] {
        auto bad = p; bad.cell_volumetric_capacity = {0.0};
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto bad = p; bad.volumetric_capacity = [](double, double, double) { return 1.0; };
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto bad = p; bad.cell_volumetric_capacity.clear();
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto bad = p; bad.cell_volumetric_capacity = {1.0, 1.0};
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto bad = p; bad.initial_temperature[0] = std::numeric_limits<double>::quiet_NaN();
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto bad = p; bad.diffusion.cell_source = {std::numeric_limits<double>::infinity()};
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto bad = p;
        bad.diffusion.dirichlet_faces = {{{0,1,2}, [](double, double, double) { return 0.0; }}};
        bad.diffusion.neumann_faces = {{{0,1,2}, [](double, double, double) { return 1.0; }}};
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto bad = p;
        bad.diffusion.robin_faces = {{{0,1,2},
            [](double, double, double) { return -1.0; },
            [](double, double, double) { return 300.0; }}};
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto bad = p;
        bad.diffusion.dirichlet_constraints = {{0,0},{0,1}};
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto options = TransientHeat3DOptions{};
        options.time_step = -0.1;
        (void)solve_transient_heat_3d(p, options);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto options = TransientHeat3DOptions{};
        options.steps = 0;
        (void)solve_transient_heat_3d(p, options);
    }));
}

NOVA_TEST(transient_heat_3d_streaming_reports_exact_states_and_retains_metadata) {
    auto p = unit_heat();
    p.cell_volumetric_capacity = {4.0};
    p.diffusion.source = [](double, double, double) { return 8.0; };
    p.initial_temperature = {300, 300, 300, 300};
    TransientHeat3DOptions opts;
    opts.start_time = 2.0;
    opts.time_step = 0.25;
    opts.steps = 4;
    opts.scheme = TransientHeat3DScheme::CrankNicolson;
    opts.retain_step_solutions = false;
    opts.solver_options.relative_tolerance = 1e-13;
    opts.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    std::vector<double> times;
    std::vector<std::vector<double>> frames;
    opts.step_observer = [&](double t, std::span<const double> values) {
        times.push_back(t);
        frames.emplace_back(values.begin(), values.end());
    };
    const auto result = solve_transient_heat_3d(p, opts);
    NOVA_REQUIRE(times.size() == 5);
    NOVA_REQUIRE(frames.size() == 5);
    NOVA_REQUIRE(result.linear_solvers.size() == 4);
    NOVA_REQUIRE(result.diagnostics.size() == 4);
    for (std::size_t step = 0; step < 5; ++step) {
        NOVA_REQUIRE_NEAR(times[step], 2.0 + 0.25 * static_cast<double>(step), 1e-14);
        NOVA_REQUIRE(frames[step].size() == 4);
        for (const double t : frames[step]) {
            NOVA_REQUIRE_NEAR(t, 300.0 + 0.5 * static_cast<double>(step), 3e-9);
        }
        if (step != 0) {
            NOVA_REQUIRE(result.linear_solvers[step - 1].converged);
            NOVA_REQUIRE(result.linear_solvers[step - 1].solution.empty());
            NOVA_REQUIRE_NEAR(result.diagnostics[step - 1].time, times[step], 1e-14);
        }
    }
    for (std::size_t node = 0; node < 4; ++node) {
        NOVA_REQUIRE_NEAR(result.temperature[node], frames.back()[node], 1e-14);
    }
    opts.step_observer = {};
    opts.retain_step_solutions = true;
    const auto legacy = solve_transient_heat_3d(p, opts);
    for (std::size_t step = 0; step < 4; ++step) {
        NOVA_REQUIRE(legacy.linear_solvers[step].solution.size() == 4);
        for (std::size_t node = 0; node < 4; ++node) {
            NOVA_REQUIRE_NEAR(legacy.linear_solvers[step].solution[node],
                              frames[step + 1][node], 1e-12);
        }
        NOVA_REQUIRE_NEAR(legacy.diagnostics[step].balance_residual,
                          result.diagnostics[step].balance_residual, 1e-12);
    }
}

NOVA_TEST(transient_heat_3d_initial_observer_receives_projected_dirichlet_state) {
    auto p = unit_heat();
    p.diffusion.dirichlet_constraints = {{0, 250.0}};
    p.initial_temperature = {10.0, 20.0, 30.0, 40.0};
    TransientHeat3DOptions opts;
    opts.start_time = 1.5;
    opts.time_step = 0.05;
    std::vector<double> initial;
    std::size_t observations = 0;
    opts.step_observer = [&](double t, std::span<const double> values) {
        if (observations++ == 0) {
            NOVA_REQUIRE_NEAR(t, 1.5, 1e-15);
            initial.assign(values.begin(), values.end());
        }
    };
    (void)solve_transient_heat_3d(p, opts);
    NOVA_REQUIRE(observations == 2);
    NOVA_REQUIRE(initial.size() == 4);
    NOVA_REQUIRE_NEAR(initial[0], 250.0, 1e-15);
    NOVA_REQUIRE_NEAR(initial[1], 20.0, 1e-15);
    NOVA_REQUIRE_NEAR(initial[2], 30.0, 1e-15);
    NOVA_REQUIRE_NEAR(initial[3], 40.0, 1e-15);
}

NOVA_TEST(transient_heat_3d_observer_exceptions_propagate) {
    auto p = unit_heat();
    TransientHeat3DOptions opts;
    opts.step_observer = [](double, std::span<const double>) {
        throw std::runtime_error("observer failure");
    };
    bool caught = false;
    try {
        (void)solve_transient_heat_3d(p, opts);
    } catch (const std::runtime_error& e) {
        caught = std::string(e.what()) == "observer failure";
    }
    NOVA_REQUIRE(caught);
}
