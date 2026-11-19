#include "test_framework.hpp"
#include "nova/fem/transient_heat_3d.hpp"

#include <cstddef>

namespace {
using namespace nova::fem;

TransientHeat3DProblem restart_case() {
    TransientHeat3DProblem p;
    p.diffusion.mesh = {
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
        {{0, 1, 2, 3}}
    };
    p.diffusion.cell_diffusivity = {5.0};
    p.diffusion.source = [](double, double, double) { return 1.5; };
    p.diffusion.cell_source = {0.5};
    p.cell_volumetric_capacity = {2.0};
    // Project the deliberately inconsistent initial temperature to 300 K.
    p.initial_temperature = {100.0, 270.0, 285.0, 290.0};
    p.diffusion.dirichlet_constraints = {{0, 300.0}};
    p.diffusion.neumann_faces.push_back({
        {0, 1, 2}, [](double, double, double) { return 0.2; }
    });
    p.diffusion.robin_faces.push_back({
        {1, 2, 3},
        [](double, double, double) { return 2.0; },
        [](double, double, double) { return 295.0; }
    });
    return p;
}
} // namespace

NOVA_TEST(transient_heat_3d_restart_matches_uninterrupted_implicit_integration) {
    // Discrete restart equivalence, not a continuum-accuracy claim.
    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                              TransientHeat3DScheme::CrankNicolson}) {
        const auto p = restart_case();
        TransientHeat3DOptions options;
        options.scheme = scheme;
        options.start_time = 4.25;
        options.time_step = 0.05;
        options.steps = 8;
        options.solver_options.preconditioner =
            nova::linalg::PreconditionerType::Jacobi;
        options.solver_options.relative_tolerance = 1e-12;

        const auto uninterrupted = solve_transient_heat_3d(p, options);
        options.steps = 3;
        const auto first = solve_transient_heat_3d(p, options);
        NOVA_REQUIRE_NEAR(first.temperature[0], 300.0, 1e-12);
        NOVA_REQUIRE_NEAR(first.final_time, 4.4, 1e-12);

        auto restarted = p;
        restarted.initial_temperature = first.temperature;
        options.start_time = first.final_time;
        options.steps = 5;
        const auto second = solve_transient_heat_3d(restarted, options);

        NOVA_REQUIRE(uninterrupted.linear_solvers.size() == 8U);
        NOVA_REQUIRE(uninterrupted.diagnostics.size() == 8U);
        NOVA_REQUIRE(first.diagnostics.size() == 3U);
        NOVA_REQUIRE(second.diagnostics.size() == 5U);
        NOVA_REQUIRE_NEAR(second.final_time, uninterrupted.final_time, 1e-12);
        NOVA_REQUIRE_NEAR(uninterrupted.final_time, 4.65, 1e-12);
        for (std::size_t node = 0; node < second.temperature.size(); ++node) {
            NOVA_REQUIRE_NEAR(
                first.temperature[node],
                uninterrupted.linear_solvers[2].solution[node], 2e-9);
            NOVA_REQUIRE_NEAR(
                second.temperature[node], uninterrupted.temperature[node], 2e-9);
        }
        for (std::size_t i = 0; i < second.diagnostics.size(); ++i) {
            const auto& resumed = second.diagnostics[i];
            const auto& continuous = uninterrupted.diagnostics[i + 3U];
            NOVA_REQUIRE_NEAR(resumed.time, continuous.time, 1e-12);
            NOVA_REQUIRE_NEAR(resumed.heat_content, continuous.heat_content, 2e-9);
            NOVA_REQUIRE_NEAR(resumed.source_energy, continuous.source_energy, 1e-11);
            NOVA_REQUIRE_NEAR(resumed.neumann_energy, continuous.neumann_energy, 1e-11);
            NOVA_REQUIRE_NEAR(resumed.robin_energy, continuous.robin_energy, 2e-9);
            NOVA_REQUIRE_NEAR(resumed.dirichlet_energy, continuous.dirichlet_energy, 2e-9);
            NOVA_REQUIRE_NEAR(resumed.balance_residual, 0.0, 1e-9);
            NOVA_REQUIRE_NEAR(continuous.balance_residual, 0.0, 1e-9);
        }
    }
}
