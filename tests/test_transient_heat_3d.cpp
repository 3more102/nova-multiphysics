#include "test_framework.hpp"
#include "nova/fem/transient_heat_3d.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

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


NOVA_TEST(transient_heat_3d_temporal_refinement_matches_be_and_cn_orders) {
    // Three fixed unit-tetra vertices leave one scalar semidiscrete mode:
    // M33=1/60, K33=1/6 => dT/dt = -10*T. Its exact solution at t=0.1
    // is 10*exp(-1), independent of either time integrator.
    auto p = unit_heat();
    p.initial_temperature = {0,0,0,10};
    p.diffusion.dirichlet_constraints = {{0,0}, {1,0}, {2,0}};
    const double exact = 10.0 * std::exp(-1.0);
    const auto error = [&](TransientHeat3DScheme scheme, std::size_t steps) {
        TransientHeat3DOptions opts;
        opts.scheme = scheme;
        opts.time_step = 0.1 / static_cast<double>(steps);
        opts.steps = steps;
        opts.solver_options.relative_tolerance = 1e-13;
        opts.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
        const auto r = solve_transient_heat_3d(p, opts);
        NOVA_REQUIRE_NEAR(r.final_time, 0.1, 1e-14);
        NOVA_REQUIRE_NEAR(r.diagnostics.back().balance_residual, 0.0, 1e-12);
        return std::abs(r.temperature[3] - exact);
    };

    const double be4 = error(TransientHeat3DScheme::BackwardEuler, 4);
    const double be8 = error(TransientHeat3DScheme::BackwardEuler, 8);
    const double cn4 = error(TransientHeat3DScheme::CrankNicolson, 4);
    const double cn8 = error(TransientHeat3DScheme::CrankNicolson, 8);
    NOVA_REQUIRE(be8 > 0.0);
    NOVA_REQUIRE(cn8 > 0.0);
    // First-order global BE error ~dt; second-order global CN error ~dt^2.
    NOVA_REQUIRE(be4 / be8 > 1.75);
    NOVA_REQUIRE(be4 / be8 < 2.1);
    NOVA_REQUIRE(cn4 / cn8 > 3.5);
    NOVA_REQUIRE(cn4 / cn8 < 4.5);
}

NOVA_TEST(transient_heat_3d_mainline_cell_tensor_matches_exact_decay) {
    // Verify the newer mainline per-tetrahedron SPD tensor path is accepted
    // by the transient assembler, not just the original spatial tensor path.
    auto p = unit_heat();
    p.diffusion.cell_diffusivity.clear();
    p.diffusion.cell_diffusivity_tensors = {
        SymmetricTensor3D{2.0,3.0,4.0,0.1,0.2,-0.1}
    };
    p.cell_volumetric_capacity = {2.0};
    p.initial_temperature = {0,0,0,10};
    p.diffusion.dirichlet_constraints = {{0,0}, {1,0}, {2,0}};
    TransientHeat3DOptions opts;
    opts.time_step = 0.05;
    opts.solver_options.relative_tolerance = 1e-13;
    opts.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    // dT/dt = -(10 * Kzz / C)*T = -20*T; one BE step halves the free DOF.
    const auto result = solve_transient_heat_3d(p, opts);
    NOVA_REQUIRE_NEAR(result.temperature[3], 5.0, 2e-12);
    NOVA_REQUIRE_NEAR(result.diagnostics[0].balance_residual, 0.0, 2e-12);
}

NOVA_TEST(transient_heat_3d_step_callback_streams_without_full_solution_history) {
    auto p = unit_heat();
    p.cell_volumetric_capacity = {4.0};
    p.diffusion.source = [](double, double, double) { return 8.0; };
    p.initial_temperature = {300.0, 300.0, 300.0, 300.0};

    TransientHeat3DOptions opts;
    opts.time_step = 0.25;
    opts.start_time = 2.0;
    opts.steps = 4;
    opts.retain_step_solution_vectors = false;

    std::size_t observed_steps = 0;
    opts.on_step = [&](std::size_t step, const std::vector<double>& temperature,
                       const TransientHeatStepDiagnostics3D& diagnostics) {
        NOVA_REQUIRE(step == observed_steps + 1);
        NOVA_REQUIRE(temperature.size() == p.diffusion.mesh.nodes.size());
        NOVA_REQUIRE_NEAR(diagnostics.time, 2.0 + 0.25 * static_cast<double>(step), 1e-14);
        NOVA_REQUIRE_NEAR(diagnostics.balance_residual, 0.0, 1e-10);
        for (const double value : temperature) {
            NOVA_REQUIRE_NEAR(value, 300.0 + 0.5 * static_cast<double>(step), 2e-9);
        }
        ++observed_steps;
    };

    const auto output = solve_transient_heat_3d(p, opts);
    NOVA_REQUIRE(observed_steps == opts.steps);
    NOVA_REQUIRE(output.linear_solvers.size() == opts.steps);
    NOVA_REQUIRE(output.diagnostics.size() == opts.steps);
    NOVA_REQUIRE_NEAR(output.final_time, 3.0, 1e-14);
    for (const auto& solver : output.linear_solvers) {
        NOVA_REQUIRE(solver.converged);
        NOVA_REQUIRE(solver.solution.empty());
    }
    for (const double value : output.temperature) {
        NOVA_REQUIRE_NEAR(value, 302.0, 2e-9);
    }

    // The default behavior remains backward-compatible: store full states.
    opts.retain_step_solution_vectors = true;
    opts.on_step = {};
    opts.steps = 1;
    const auto default_output = solve_transient_heat_3d(p, opts);
    NOVA_REQUIRE(default_output.linear_solvers.size() == 1);
    NOVA_REQUIRE(default_output.linear_solvers[0].solution.size() ==
                 p.diffusion.mesh.nodes.size());
    for (std::size_t i = 0; i < p.diffusion.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(default_output.linear_solvers[0].solution[i],
                          default_output.temperature[i], 1e-14);
    }
}
