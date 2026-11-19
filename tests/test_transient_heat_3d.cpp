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

NOVA_TEST(transient_heat_3d_integrated_cell_tensors_preserve_source_energy) {
    auto p = unit_heat();
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 1, 1, 1);
    p.diffusion.mesh = box.mesh;
    p.diffusion.cell_diffusivity.clear();
    p.diffusion.cell_diffusivity_tensors.assign(box.mesh.elements.size(),
        SymmetricTensor3D{3.0, 2.0, 1.0, 0.2, 0.1, 0.0});
    for (std::size_t i = 1; i < p.diffusion.cell_diffusivity_tensors.size(); i += 2) {
        p.diffusion.cell_diffusivity_tensors[i] =
            SymmetricTensor3D{4.0, 2.0, 1.5, -0.2, 0.0, 0.1};
    }
    p.diffusion.source = [](double, double, double) { return 8.0; };
    p.cell_volumetric_capacity.assign(box.mesh.elements.size(), 4.0);
    p.initial_temperature.assign(box.mesh.nodes.size(), 300.0);
    TransientHeat3DOptions options;
    options.time_step = 0.25;
    options.steps = 6;
    options.scheme = TransientHeat3DScheme::CrankNicolson;
    options.solver_options.relative_tolerance = 1e-12;
    options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    const auto result = solve_transient_heat_3d(p, options);
    NOVA_REQUIRE(result.diagnostics.size() == 6U);
    for (double temperature : result.temperature) {
        NOVA_REQUIRE_NEAR(temperature, 303.0, 2e-8);
    }
    for (const auto& d : result.diagnostics) {
        NOVA_REQUIRE_NEAR(d.source_energy, 2.0, 1e-10);
        NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1e-8);
    }
    NOVA_REQUIRE_NEAR(result.diagnostics.back().heat_content, 1212.0, 2e-7);
}

NOVA_TEST(transient_heat_3d_integrated_cell_tensors_validate_input) {
    auto p = unit_heat();
    p.diffusion.cell_diffusivity.clear();
    p.diffusion.cell_diffusivity_tensors = {
        SymmetricTensor3D{2.0, 3.0, 4.0, 0.0, 0.0, 0.0}
    };
    const auto valid = solve_transient_heat_3d(p);
    NOVA_REQUIRE(valid.temperature.size() == 4U);
    NOVA_REQUIRE(rejects([&] {
        auto bad = p;
        bad.diffusion.cell_diffusivity_tensors[0].xx = -1.0;
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto bad = p;
        bad.diffusion.cell_diffusivity_tensors.push_back(
            SymmetricTensor3D{2.0, 3.0, 4.0, 0.0, 0.0, 0.0});
        (void)solve_transient_heat_3d(bad);
    }));
    NOVA_REQUIRE(rejects([&] {
        auto bad = p;
        bad.diffusion.cell_diffusivity = {1.0};
        (void)solve_transient_heat_3d(bad);
    }));
}


NOVA_TEST(transient_heat_3d_temporal_convergence_matches_exact_decay) {
    // A single unit tetrahedron with only vertex 3 free has
    // M_33 = 1/60 and K_33 = 1/6, hence dT_3/dt = -10 T_3.
    // The exact final temperature at t=0.2 with T_3(0)=10 is
    // 10*exp(-2). Halving dt should reduce BE error by ~2 and
    // CN error by ~4. This separately verifies temporal order,
    // not just one-step algebra or global heat conservation.
    const double exact = 10.0 * std::exp(-2.0);
    const auto temporal_error = [&](TransientHeat3DScheme scheme,
                                    double dt, std::size_t steps) {
        auto p = unit_heat();
        p.initial_temperature = {0.0, 0.0, 0.0, 10.0};
        p.diffusion.dirichlet_constraints = {{0,0}, {1,0}, {2,0}};
        TransientHeat3DOptions options;
        options.scheme = scheme;
        options.start_time = 3.0;
        options.time_step = dt;
        options.steps = steps;
        options.solver_options.relative_tolerance = 1e-13;
        options.solver_options.preconditioner =
            nova::linalg::PreconditionerType::Jacobi;
        const auto result = solve_transient_heat_3d(p, options);
        NOVA_REQUIRE_NEAR(result.final_time, 3.2, 1e-14);
        NOVA_REQUIRE(result.diagnostics.size() == steps);
        NOVA_REQUIRE_NEAR(result.diagnostics.back().time, 3.2, 1e-14);
        for (const auto& d : result.diagnostics) {
            NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1e-11);
        }
        return std::abs(result.temperature[3] - exact);
    };

    const double be_coarse = temporal_error(
        TransientHeat3DScheme::BackwardEuler, 0.05, 4);
    const double be_fine = temporal_error(
        TransientHeat3DScheme::BackwardEuler, 0.025, 8);
    NOVA_REQUIRE(be_fine > 0.0);
    NOVA_REQUIRE(be_coarse / be_fine > 1.75);
    NOVA_REQUIRE(be_coarse / be_fine < 2.2);

    const double cn_coarse = temporal_error(
        TransientHeat3DScheme::CrankNicolson, 0.05, 4);
    const double cn_fine = temporal_error(
        TransientHeat3DScheme::CrankNicolson, 0.025, 8);
    NOVA_REQUIRE(cn_fine > 0.0);
    NOVA_REQUIRE(cn_coarse / cn_fine > 3.5);
    NOVA_REQUIRE(cn_coarse / cn_fine < 4.5);
    NOVA_REQUIRE(cn_fine < be_fine);
}


NOVA_TEST(transient_heat_3d_time_integration_exhibits_first_and_second_order) {
    // A single tetra with three vertices held at zero leaves one discrete
    // degree of freedom. Consistent P1 mass gives M_33 = 1/60, while
    // conductivity gives K_33 = 1/6, hence dT_3/dt = -10 T_3.
    // At t = 0.2 s the exact semidiscrete value is exp(-2).
    auto problem = unit_heat();
    problem.initial_temperature = {0.0, 0.0, 0.0, 1.0};
    problem.diffusion.dirichlet_constraints = {{0, 0.0}, {1, 0.0}, {2, 0.0}};
    const double exact = std::exp(-2.0);

    const auto error_at = [&](TransientHeat3DScheme scheme, double dt,
                              std::size_t steps) {
        TransientHeat3DOptions options;
        options.scheme = scheme;
        options.time_step = dt;
        options.steps = steps;
        options.solver_options.relative_tolerance = 1e-13;
        options.solver_options.preconditioner =
            nova::linalg::PreconditionerType::Jacobi;
        const auto result = solve_transient_heat_3d(problem, options);
        NOVA_REQUIRE_NEAR(result.final_time, 0.2, 1e-14);
        NOVA_REQUIRE_NEAR(result.diagnostics.back().balance_residual,
                          0.0, 1e-11);
        return std::abs(result.temperature[3] - exact);
    };

    const double be_coarse = error_at(TransientHeat3DScheme::BackwardEuler,
                                      0.05, 4);
    const double be_fine = error_at(TransientHeat3DScheme::BackwardEuler,
                                    0.025, 8);
    const double cn_coarse = error_at(TransientHeat3DScheme::CrankNicolson,
                                      0.05, 4);
    const double cn_fine = error_at(TransientHeat3DScheme::CrankNicolson,
                                    0.025, 8);

    // Halving dt should reduce BE error by approximately 2x and CN by 4x.
    // Conservative gates tolerate floating-point solver perturbations.
    NOVA_REQUIRE(be_coarse > 1.8 * be_fine);
    NOVA_REQUIRE(cn_coarse > 3.5 * cn_fine);
    NOVA_REQUIRE(cn_fine < be_fine);
}

NOVA_TEST(transient_heat_3d_streams_projected_initial_and_each_solved_step) {
    auto p = unit_heat();
    p.cell_volumetric_capacity = {4.0};
    p.diffusion.source = [](double, double, double) { return 8.0; };
    p.initial_temperature = {123, 300, 300, 300};
    p.diffusion.dirichlet_constraints = {{0, 300.0}};
    TransientHeat3DOptions options;
    options.steps = 3;
    options.time_step = 0.25;
    options.start_time = 2.0;
    options.retain_step_solutions = false;
    options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    std::vector<double> times;
    std::vector<std::vector<double>> snapshots;
    options.step_observer = [&](double t, const std::vector<double>& values) {
        times.push_back(t);
        snapshots.push_back(values);
    };
    const auto result = solve_transient_heat_3d(p, options);
    NOVA_REQUIRE(times.size() == options.steps + 1);
    NOVA_REQUIRE(snapshots.size() == options.steps + 1);
    NOVA_REQUIRE_NEAR(times.front(), options.start_time, 1e-14);
    NOVA_REQUIRE_NEAR(times.back(), result.final_time, 1e-14);
    NOVA_REQUIRE_NEAR(snapshots.front()[0], 300.0, 1e-14);
    NOVA_REQUIRE(result.linear_solvers.size() == options.steps);
    NOVA_REQUIRE(result.diagnostics.size() == options.steps);
    for (const auto& solver : result.linear_solvers) {
        NOVA_REQUIRE(solver.solution.empty());
        NOVA_REQUIRE(solver.converged);
    }
    for (std::size_t i = 0; i < result.temperature.size(); ++i) {
        NOVA_REQUIRE_NEAR(snapshots.back()[i], result.temperature[i], 1e-13);
    }
}

NOVA_TEST(transient_heat_3d_observer_errors_propagate) {
    auto p = unit_heat();
    TransientHeat3DOptions options;
    options.steps = 3;
    std::size_t count = 0;
    options.step_observer = [&](double, const std::vector<double>&) {
        if (++count == 2) throw std::domain_error("Output failed.");
    };
    bool propagated = false;
    try {
        (void)solve_transient_heat_3d(p, options);
    } catch (const std::domain_error&) {
        propagated = true;
    }
    NOVA_REQUIRE(propagated);
    NOVA_REQUIRE(count == 2);
}

NOVA_TEST(transient_heat_3d_rejects_nonadvancing_step_time) {
    auto p = unit_heat();
    TransientHeat3DOptions options;
    options.start_time = 1.0e20;
    options.time_step = 1.0e3;
    options.steps = 100;
    NOVA_REQUIRE(rejects([&] { (void)solve_transient_heat_3d(p, options); }));
}
