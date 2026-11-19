#include "test_framework.hpp"
#include "nova/fem/transient_heat_3d.hpp"

#include <cstddef>
#include <vector>

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
    // The solver projects this deliberately inconsistent field to 300 K at node 0.
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

void capture_frame(
    std::vector<std::vector<double>>& frames,
    const std::vector<double>& temperature
) {
    frames.push_back(temperature);
}
} // namespace

NOVA_TEST(transient_heat_3d_restart_matches_uninterrupted_implicit_integration) {
    // This verifies discrete restart equivalence, not continuum accuracy.
    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                              TransientHeat3DScheme::CrankNicolson}) {
        const auto p = restart_case();
        TransientHeat3DOptions options;
        options.scheme = scheme;
        options.start_time = 4.25;
        options.time_step = 0.05;
        options.steps = 8;
        options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
        options.solver_options.relative_tolerance = 1e-12;

        std::vector<std::vector<double>> uninterrupted_frames;
        options.on_step = [&](std::size_t, double, const std::vector<double>& field) {
            capture_frame(uninterrupted_frames, field);
        };
        const auto uninterrupted = solve_transient_heat_3d(p, options);

        std::vector<std::vector<double>> first_frames;
        options.steps = 3;
        options.on_step = [&](std::size_t, double, const std::vector<double>& field) {
            capture_frame(first_frames, field);
        };
        const auto first = solve_transient_heat_3d(p, options);
        NOVA_REQUIRE_NEAR(first.temperature[0], 300.0, 1e-12);
        NOVA_REQUIRE_NEAR(first.final_time, 4.4, 1e-12);

        auto restarted = p;
        restarted.initial_temperature = first.temperature;
        options.start_time = first.final_time;
        options.steps = 5;
        std::vector<std::vector<double>> resumed_frames;
        options.on_step = [&](std::size_t, double, const std::vector<double>& field) {
            capture_frame(resumed_frames, field);
        };
        const auto second = solve_transient_heat_3d(restarted, options);

        NOVA_REQUIRE(uninterrupted.linear_solvers.size() == 8U);
        NOVA_REQUIRE(uninterrupted.diagnostics.size() == 8U);
        NOVA_REQUIRE(first.diagnostics.size() == 3U);
        NOVA_REQUIRE(second.diagnostics.size() == 5U);
        NOVA_REQUIRE(uninterrupted_frames.size() == 9U);
        NOVA_REQUIRE(first_frames.size() == 4U);
        NOVA_REQUIRE(resumed_frames.size() == 6U);
        NOVA_REQUIRE_NEAR(second.final_time, uninterrupted.final_time, 1e-12);
        NOVA_REQUIRE_NEAR(uninterrupted.final_time, 4.65, 1e-12);

        for (std::size_t node = 0; node < second.temperature.size(); ++node) {
            NOVA_REQUIRE_NEAR(first.temperature[node], uninterrupted_frames[3][node], 2e-9);
            NOVA_REQUIRE_NEAR(first.temperature[node], resumed_frames[0][node], 2e-9);
            NOVA_REQUIRE_NEAR(second.temperature[node], uninterrupted.temperature[node], 2e-9);
            NOVA_REQUIRE_NEAR(second.temperature[node], uninterrupted_frames[8][node], 2e-9);
        }
        for (std::size_t step = 0; step < 5U; ++step) {
            for (std::size_t node = 0; node < second.temperature.size(); ++node) {
                NOVA_REQUIRE_NEAR(
                    resumed_frames[step + 1U][node],
                    uninterrupted_frames[step + 4U][node], 2e-9);
            }
            const auto& resumed = second.diagnostics[step];
            const auto& continuous = uninterrupted.diagnostics[step + 3U];
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
