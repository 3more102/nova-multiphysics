#include "test_framework.hpp"
#include "nova/fem/transient_heat_3d.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <string>
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


// Independent temporal-convergence check against the *semidiscrete* exact
// decay solution, not against the solver's own update formulas. Three fixed
// vertices leave one free DOF: M=1/60, K=1/6, hence T'(t)=-10*T(t).
// On a fixed spatial mesh BE must converge at order one and CN at order two.
NOVA_TEST(transient_heat_3d_temporal_convergence_matches_exact_decay) {
    auto p = unit_heat();
    p.initial_temperature = {0.0, 0.0, 0.0, 10.0};
    p.diffusion.dirichlet_constraints = {{0,0}, {1,0}, {2,0}};
    constexpr double finish_time = 0.2;
    const double exact = 10.0 * std::exp(-10.0 * finish_time);
    const std::array<std::size_t, 3> step_counts{10U, 20U, 40U};

    TransientHeat3DOptions options;
    options.solver_options.relative_tolerance = 1e-13;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::Jacobi;
    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                               TransientHeat3DScheme::CrankNicolson}) {
        options.scheme = scheme;
        std::array<double, 3> errors{};
        for (std::size_t i = 0; i < step_counts.size(); ++i) {
            options.steps = step_counts[i];
            options.time_step = finish_time / static_cast<double>(options.steps);
            const auto result = solve_transient_heat_3d(p, options);
            NOVA_REQUIRE_NEAR(result.final_time, finish_time, 1e-14);
            errors[i] = std::abs(result.temperature[3] - exact);
            NOVA_REQUIRE(std::isfinite(errors[i]) && errors[i] > 0.0);
        }
        NOVA_REQUIRE(errors[0] > errors[1] && errors[1] > errors[2]);
        for (std::size_t i = 0; i + 1 < errors.size(); ++i) {
            const double reduction = errors[i] / errors[i + 1];
            if (scheme == TransientHeat3DScheme::BackwardEuler) {
                NOVA_REQUIRE(reduction > 1.8 && reduction < 2.2);
            } else {
                NOVA_REQUIRE(reduction > 3.8 && reduction < 4.2);
            }
        }
    }
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


NOVA_TEST(transient_heat_3d_streaming_mode_matches_full_history_for_both_schemes) {
    auto p = unit_heat();
    p.cell_volumetric_capacity = {4.0};
    p.diffusion.source = [](double, double, double) { return 8.0; };
    p.initial_temperature = {300.0, 300.0, 300.0, 300.0};

    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                              TransientHeat3DScheme::CrankNicolson}) {
        TransientHeat3DOptions options;
        options.scheme = scheme;
        options.start_time = 2.0;
        options.time_step = 0.125;
        options.steps = 10;
        options.store_step_history = false;
        options.solver_options.relative_tolerance = 1e-13;
        options.solver_options.preconditioner =
            nova::linalg::PreconditionerType::Jacobi;
        std::size_t callback_count = 0;
        std::vector<double> last_callback_field;

        options.on_step = [&](std::size_t step, double time,
                              std::span<const double> temperatures,
                              const TransientHeatStepDiagnostics3D& diag) {
            ++callback_count;
            NOVA_REQUIRE(step == callback_count);
            NOVA_REQUIRE_NEAR(time, 2.0 + step * 0.125, 1e-14);
            NOVA_REQUIRE_NEAR(diag.time, time, 1e-14);
            NOVA_REQUIRE_NEAR(diag.source_energy, 1.0 / 6.0, 1e-13);
            NOVA_REQUIRE_NEAR(diag.balance_residual, 0.0, 1e-10);
            NOVA_REQUIRE(temperatures.size() == 4U);
            for (double value : temperatures) {
                NOVA_REQUIRE_NEAR(value, 300.0 + 2.0 * step * 0.125, 2e-9);
            }
            last_callback_field.assign(temperatures.begin(), temperatures.end());
        };

        const auto streamed = solve_transient_heat_3d(p, options);
        NOVA_REQUIRE(callback_count == options.steps);
        NOVA_REQUIRE(streamed.linear_solvers.empty());
        NOVA_REQUIRE(streamed.diagnostics.empty());
        NOVA_REQUIRE_NEAR(streamed.final_time, 3.25, 1e-14);
        NOVA_REQUIRE(streamed.temperature == last_callback_field);

        options.store_step_history = true;
        options.on_step = {};
        const auto retained = solve_transient_heat_3d(p, options);
        NOVA_REQUIRE(retained.linear_solvers.size() == options.steps);
        NOVA_REQUIRE(retained.diagnostics.size() == options.steps);
        for (std::size_t i = 0; i < streamed.temperature.size(); ++i) {
            NOVA_REQUIRE_NEAR(streamed.temperature[i], retained.temperature[i], 1e-13);
            NOVA_REQUIRE_NEAR(retained.linear_solvers.back().solution[i],
                              retained.temperature[i], 1e-13);
        }
    }
}

NOVA_TEST(transient_heat_3d_stream_callback_exception_aborts_the_solve) {
    auto p = unit_heat();
    std::size_t callback_count = 0;
    TransientHeat3DOptions options;
    options.steps = 4;
    options.store_step_history = false;
    options.on_step = [&](std::size_t step, double,
                          std::span<const double>,
                          const TransientHeatStepDiagnostics3D&) {
        ++callback_count;
        if (step == 2U) throw std::runtime_error("stop streaming");
    };
    bool caught = false;
    try {
        (void)solve_transient_heat_3d(p, options);
    } catch (const std::runtime_error& error) {
        caught = std::string(error.what()) == "stop streaming";
    }
    NOVA_REQUIRE(caught);
    NOVA_REQUIRE(callback_count == 2U);
}
