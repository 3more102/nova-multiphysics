#include "test_framework.hpp"
#include "nova/fem/adaptive_heat_3d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace nova::fem;

TransientHeat3DProblem tetra_problem() {
    TransientHeat3DProblem problem;
    problem.diffusion.mesh = {
        {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}
    };
    problem.diffusion.cell_diffusivity = {1.0};
    problem.cell_volumetric_capacity = {1.0};
    problem.initial_temperature = {0.0, 0.0, 0.0, 10.0};
    problem.diffusion.dirichlet_constraints = {{0,0}, {1,0}, {2,0}};
    return problem;
}

template<class F>
bool rejects_invalid(F&& function) {
    try { function(); }
    catch (const std::invalid_argument&) { return true; }
    return false;
}

} // namespace

NOVA_TEST(adaptive_heat_3d_exact_uniform_source_and_energy) {
    auto p = tetra_problem();
    p.diffusion.dirichlet_constraints.clear();
    p.cell_volumetric_capacity = {4.0};
    p.diffusion.source = [](double, double, double) { return 8.0; };
    p.initial_temperature = {300.0,300.0,300.0,300.0};
    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                              TransientHeat3DScheme::CrankNicolson}) {
        AdaptiveHeat3DOptions options;
        options.scheme = scheme;
        options.start_time = 2.5;
        options.end_time = 3.5;
        options.initial_step = 0.25;
        options.max_step = 0.4;
        options.absolute_tolerance = 1e-7;
        options.relative_tolerance = 0.0;
        options.solver_options.relative_tolerance = 1e-12;
        const auto result = solve_adaptive_heat_3d(p, options);
        NOVA_REQUIRE_NEAR(result.final_time, 3.5, 1e-13);
        NOVA_REQUIRE(result.accepted_steps >= 3U);
        NOVA_REQUIRE(result.rejected_steps == 0U);
        NOVA_REQUIRE(result.diagnostics.size() == result.accepted_steps);
        NOVA_REQUIRE(result.substep_solvers.size() == 2U * result.accepted_steps);
        NOVA_REQUIRE(result.attempts.size() == result.accepted_steps);
        for (const double temperature : result.temperature) {
            NOVA_REQUIRE_NEAR(temperature, 302.0, 3e-8);
        }
        double total_source = 0.0;
        double previous_time = options.start_time;
        for (const auto& d : result.diagnostics) {
            NOVA_REQUIRE(d.time > previous_time);
            previous_time = d.time;
            total_source += d.source_energy;
            NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 2e-9);
        }
        NOVA_REQUIRE_NEAR(total_source, 8.0 / 6.0, 2e-11);
        NOVA_REQUIRE_NEAR(result.diagnostics.back().heat_content,
                          4.0 * 302.0 / 6.0, 3e-8);
    }
}


NOVA_TEST(adaptive_heat_3d_time_dependent_uniform_source_uses_cn_quadrature) {
    auto problem = tetra_problem();
    problem.diffusion.dirichlet_constraints.clear();
    problem.diffusion.source = {};
    problem.diffusion.cell_source.clear();
    problem.cell_volumetric_capacity = {2.0};
    problem.initial_temperature.assign(4U, 300.0);
    problem.time_dependent_source =
        [](double, double, double, double time) { return time; };

    AdaptiveHeat3DOptions options;
    options.scheme = TransientHeat3DScheme::CrankNicolson;
    options.start_time = 2.0;
    options.end_time = 3.0;
    options.initial_step = 0.1;
    options.min_step = 0.1;
    options.max_step = 0.1;
    options.absolute_tolerance = 1.0e-6;
    options.relative_tolerance = 0.0;
    options.solver_options.relative_tolerance = 1.0e-12;

    const auto result = solve_adaptive_heat_3d(problem, options);
    NOVA_REQUIRE(result.accepted_steps == 10U);
    NOVA_REQUIRE(result.rejected_steps == 0U);
    NOVA_REQUIRE_NEAR(result.final_time, 3.0, 1.0e-14);
    for (const double temperature : result.temperature) {
        NOVA_REQUIRE_NEAR(temperature, 301.25, 2.0e-9);
    }
    double total_source_energy = 0.0;
    for (const auto& diagnostic : result.diagnostics) {
        total_source_energy += diagnostic.source_energy;
        NOVA_REQUIRE_NEAR(diagnostic.balance_residual, 0.0, 2.0e-9);
    }
    NOVA_REQUIRE_NEAR(total_source_energy, 5.0 / 12.0, 2.0e-11);
}

NOVA_TEST(adaptive_heat_3d_tolerance_controls_decay_error) {
    const auto p = tetra_problem();
    AdaptiveHeat3DOptions loose;
    loose.scheme = TransientHeat3DScheme::BackwardEuler;
    loose.end_time = 0.2;
    loose.initial_step = 0.2;
    loose.max_step = 0.2;
    loose.min_step = 1e-5;
    loose.absolute_tolerance = 0.15;
    loose.relative_tolerance = 0.0;
    loose.solver_options.relative_tolerance = 1e-13;
    auto tight = loose;
    tight.absolute_tolerance = 0.002;

    const auto coarse = solve_adaptive_heat_3d(p, loose);
    const auto fine = solve_adaptive_heat_3d(p, tight);
    const double exact = 10.0 * std::exp(-2.0);
    const double coarse_error = std::abs(coarse.temperature[3] - exact);
    const double fine_error = std::abs(fine.temperature[3] - exact);
    NOVA_REQUIRE(fine.rejected_steps > 0U);
    NOVA_REQUIRE(fine.accepted_steps > coarse.accepted_steps);
    NOVA_REQUIRE(fine_error < coarse_error);
    NOVA_REQUIRE(fine_error < 0.05);
    NOVA_REQUIRE_NEAR(fine.final_time, 0.2, 1e-14);
    NOVA_REQUIRE_NEAR(fine.temperature[0], 0.0, 1e-14);
    for (const auto& a : fine.attempts) {
        NOVA_REQUIRE(std::isfinite(a.normalized_error));
        NOVA_REQUIRE(a.accepted == (a.normalized_error <= 1.0));
    }
}

NOVA_TEST(adaptive_heat_3d_cn_restart_and_conservation) {
    auto p = tetra_problem();
    p.diffusion.dirichlet_constraints = {{0,300.0}, {1,300.0}, {2,300.0}};
    p.initial_temperature = {100.0, 300.0, 300.0, 310.0};
    p.diffusion.neumann_faces.push_back({
        {0,1,2}, [](double, double, double) { return 0.1; }
    });
    AdaptiveHeat3DOptions options;
    options.scheme = TransientHeat3DScheme::CrankNicolson;
    options.start_time = 1.0;
    options.end_time = 1.2;
    options.initial_step = 0.1;
    options.max_step = 0.1;
    options.absolute_tolerance = 0.01;
    options.relative_tolerance = 0.0;
    options.solver_options.relative_tolerance = 1e-12;
    const auto single = solve_adaptive_heat_3d(p, options);
    options.end_time = 1.1;
    const auto first = solve_adaptive_heat_3d(p, options);
    p.initial_temperature = first.temperature;
    options.start_time = first.final_time;
    options.end_time = 1.2;
    const auto second = solve_adaptive_heat_3d(p, options);
    NOVA_REQUIRE_NEAR(single.final_time, second.final_time, 1e-14);
    NOVA_REQUIRE_NEAR(second.temperature[0], 300.0, 1e-13);
    NOVA_REQUIRE(std::abs(single.temperature[3] - second.temperature[3]) < 0.03);
    for (const auto& d : single.diagnostics) {
        NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1e-9);
    }
    for (const auto& d : second.diagnostics) {
        NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1e-9);
    }
}

NOVA_TEST(adaptive_heat_3d_rejects_invalid_controls_and_attempt_limit) {
    auto p = tetra_problem();
    const AdaptiveHeat3DOptions defaults;
    NOVA_REQUIRE(rejects_invalid([&] {
        auto o = defaults; o.end_time = o.start_time;
        (void)solve_adaptive_heat_3d(p,o);
    }));
    NOVA_REQUIRE(rejects_invalid([&] {
        auto o = defaults; o.max_step = 1e-9;
        (void)solve_adaptive_heat_3d(p,o);
    }));
    NOVA_REQUIRE(rejects_invalid([&] {
        auto o = defaults; o.absolute_tolerance = 0;
        (void)solve_adaptive_heat_3d(p,o);
    }));
    NOVA_REQUIRE(rejects_invalid([&] {
        auto o = defaults; o.relative_tolerance = -1;
        (void)solve_adaptive_heat_3d(p,o);
    }));
    NOVA_REQUIRE(rejects_invalid([&] {
        auto o = defaults; o.initial_step = std::numeric_limits<double>::infinity();
        (void)solve_adaptive_heat_3d(p,o);
    }));
    NOVA_REQUIRE(rejects_invalid([&] {
        auto o = defaults; o.max_attempts = 0;
        (void)solve_adaptive_heat_3d(p,o);
    }));
    NOVA_REQUIRE(rejects_invalid([&] {
        auto o = defaults; o.scheme = static_cast<TransientHeat3DScheme>(100);
        (void)solve_adaptive_heat_3d(p,o);
    }));
    NOVA_REQUIRE(rejects_invalid([&] {
        auto q = p; q.initial_temperature.pop_back();
        (void)solve_adaptive_heat_3d(q, defaults);
    }));
    AdaptiveHeat3DOptions capped;
    capped.end_time = 0.2;
    capped.initial_step = 0.2;
    capped.max_step = 0.2;
    capped.absolute_tolerance = 1e-8;
    capped.relative_tolerance = 0.0;
    capped.max_attempts = 1U;
    bool capped_failure = false;
    try { (void)solve_adaptive_heat_3d(p, capped); }
    catch (const std::runtime_error&) { capped_failure = true; }
    NOVA_REQUIRE(capped_failure);
}


NOVA_TEST(adaptive_heat_3d_observer_emits_projected_initial_and_accepted_fields_only) {
    auto p = tetra_problem();
    p.initial_temperature = {5.0, 6.0, 7.0, 10.0};

    AdaptiveHeat3DOptions options;
    options.scheme = TransientHeat3DScheme::BackwardEuler;
    options.start_time = 3.0;
    options.end_time = 3.2;
    options.initial_step = 0.2;
    options.min_step = 1.0e-5;
    options.max_step = 0.2;
    options.absolute_tolerance = 0.002;
    options.relative_tolerance = 0.0;
    options.solver_options.relative_tolerance = 1.0e-13;

    std::vector<double> observed_times;
    std::vector<std::vector<double>> observed_fields;
    options.on_step = [&](std::size_t index, double time,
                          const std::vector<double>& temperature) {
        NOVA_REQUIRE(index == observed_times.size());
        observed_times.push_back(time);
        observed_fields.push_back(temperature);
    };

    const auto result = solve_adaptive_heat_3d(p, options);
    NOVA_REQUIRE(result.rejected_steps > 0U);
    NOVA_REQUIRE(observed_times.size() == result.accepted_steps + 1U);
    NOVA_REQUIRE(observed_fields.size() == observed_times.size());
    NOVA_REQUIRE_NEAR(observed_times.front(), options.start_time, 1.0e-14);
    NOVA_REQUIRE_NEAR(observed_fields.front()[0], 0.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(observed_fields.front()[1], 0.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(observed_fields.front()[2], 0.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(observed_fields.front()[3], 10.0, 1.0e-14);
    for (std::size_t i = 1U; i < observed_times.size(); ++i) {
        NOVA_REQUIRE(observed_times[i] > observed_times[i - 1U]);
        NOVA_REQUIRE_NEAR(observed_times[i], result.diagnostics[i - 1U].time, 1.0e-13);
    }
    NOVA_REQUIRE_NEAR(observed_times.back(), result.final_time, 1.0e-13);
    for (std::size_t i = 0U; i < result.temperature.size(); ++i) {
        NOVA_REQUIRE_NEAR(observed_fields.back()[i], result.temperature[i], 1.0e-13);
    }
}

NOVA_TEST(adaptive_heat_3d_observer_exceptions_propagate) {
    auto p = tetra_problem();
    AdaptiveHeat3DOptions options;
    options.on_step = [](std::size_t, double, const std::vector<double>&) {
        throw std::runtime_error("initial observer aborted");
    };
    bool initial_propagated = false;
    try {
        (void)solve_adaptive_heat_3d(p, options);
    } catch (const std::runtime_error& error) {
        initial_propagated = std::string(error.what()) == "initial observer aborted";
    }
    NOVA_REQUIRE(initial_propagated);

    auto uniform = tetra_problem();
    uniform.diffusion.dirichlet_constraints.clear();
    uniform.diffusion.source = [](double, double, double) { return 2.0; };
    uniform.initial_temperature = {10.0, 10.0, 10.0, 10.0};
    options.start_time = 0.0;
    options.end_time = 0.1;
    options.initial_step = 0.1;
    options.min_step = 0.1;
    options.max_step = 0.1;
    options.absolute_tolerance = 1.0e-8;
    options.relative_tolerance = 0.0;
    options.solver_options.relative_tolerance = 1.0e-13;
    options.on_step = [](std::size_t index, double, const std::vector<double>&) {
        if (index == 1U) throw std::runtime_error("accepted observer aborted");
    };
    bool accepted_propagated = false;
    try {
        (void)solve_adaptive_heat_3d(uniform, options);
    } catch (const std::runtime_error& error) {
        accepted_propagated = std::string(error.what()) == "accepted observer aborted";
    }
    NOVA_REQUIRE(accepted_propagated);
}


NOVA_TEST(adaptive_heat_3d_reuses_operators_and_ic0_factors_for_repeated_dt) {
    auto problem = tetra_problem();
    problem.diffusion.dirichlet_constraints.clear();
    problem.diffusion.source = [](double, double, double) { return 1.0; };
    problem.initial_temperature = {1.0, 1.0, 1.0, 1.0};

    AdaptiveHeat3DOptions options;
    options.scheme = TransientHeat3DScheme::CrankNicolson;
    options.start_time = 0.0;
    options.end_time = 0.5;
    options.initial_step = 0.125;
    options.min_step = 0.125;
    options.max_step = 0.125;
    options.absolute_tolerance = 1.0e-7;
    options.relative_tolerance = 0.0;
    options.solver_options.relative_tolerance = 1.0e-12;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::IncompleteCholesky0;

    const auto ic0 = solve_adaptive_heat_3d(problem, options);
    NOVA_REQUIRE(ic0.accepted_steps == 4U);
    NOVA_REQUIRE(ic0.rejected_steps == 0U);
    NOVA_REQUIRE(ic0.full_step_operator_builds == 1U);
    NOVA_REQUIRE(ic0.half_step_operator_builds == 1U);
    NOVA_REQUIRE(ic0.full_step_ic0_factor_builds == 1U);
    NOVA_REQUIRE(ic0.half_step_ic0_factor_builds == 1U);
    NOVA_REQUIRE_NEAR(ic0.final_time, 0.5, 1.0e-14);
    for (const double temperature : ic0.temperature) {
        NOVA_REQUIRE_NEAR(temperature, 1.5, 1.0e-10);
    }

    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::Jacobi;
    const auto jacobi = solve_adaptive_heat_3d(problem, options);
    NOVA_REQUIRE(jacobi.full_step_operator_builds == 1U);
    NOVA_REQUIRE(jacobi.half_step_operator_builds == 1U);
    NOVA_REQUIRE(jacobi.full_step_ic0_factor_builds == 0U);
    NOVA_REQUIRE(jacobi.half_step_ic0_factor_builds == 0U);
    NOVA_REQUIRE(jacobi.accepted_steps == ic0.accepted_steps);
    for (std::size_t i = 0; i < ic0.temperature.size(); ++i) {
        NOVA_REQUIRE_NEAR(ic0.temperature[i], jacobi.temperature[i], 1.0e-11);
    }
}

NOVA_TEST(adaptive_heat_3d_rebuilds_cached_operators_when_dt_changes) {
    auto problem = tetra_problem();
    AdaptiveHeat3DOptions options;
    options.scheme = TransientHeat3DScheme::BackwardEuler;
    options.start_time = 0.0;
    options.end_time = 0.2;
    options.initial_step = 0.2;
    options.min_step = 1.0e-5;
    options.max_step = 0.2;
    options.absolute_tolerance = 0.002;
    options.relative_tolerance = 0.0;
    options.solver_options.relative_tolerance = 1.0e-13;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::IncompleteCholesky0;

    const auto result = solve_adaptive_heat_3d(problem, options);
    NOVA_REQUIRE(result.rejected_steps > 0U);
    NOVA_REQUIRE(result.accepted_steps > 0U);
    NOVA_REQUIRE(result.full_step_operator_builds > 1U);
    NOVA_REQUIRE(result.half_step_operator_builds > 1U);
    NOVA_REQUIRE(result.full_step_ic0_factor_builds > 1U);
    NOVA_REQUIRE(result.half_step_ic0_factor_builds > 1U);
    NOVA_REQUIRE_NEAR(result.final_time, options.end_time, 1.0e-14);
}


NOVA_TEST(adaptive_heat_3d_evaluates_static_material_and_load_fields_once_per_study) {
    auto problem = tetra_problem();
    problem.diffusion.dirichlet_constraints.clear();
    problem.diffusion.cell_diffusivity.clear();
    problem.cell_volumetric_capacity.clear();
    problem.initial_temperature = {1.0, 1.0, 1.0, 1.0};

    std::size_t conductivity_samples = 0U;
    std::size_t capacity_samples = 0U;
    std::size_t source_samples = 0U;
    problem.diffusion.diffusivity =
        [&](double, double, double) {
            ++conductivity_samples;
            return 1.0;
        };
    problem.volumetric_capacity =
        [&](double, double, double) {
            ++capacity_samples;
            return 1.0;
        };
    problem.diffusion.source =
        [&](double, double, double) {
            ++source_samples;
            return 1.0;
        };

    AdaptiveHeat3DOptions options;
    options.scheme = TransientHeat3DScheme::CrankNicolson;
    options.end_time = 0.5;
    options.initial_step = 0.125;
    options.min_step = 0.125;
    options.max_step = 0.125;
    options.absolute_tolerance = 1.0e-7;
    options.relative_tolerance = 0.0;

    const auto result = solve_adaptive_heat_3d(problem, options);
    NOVA_REQUIRE(result.accepted_steps == 4U);
    NOVA_REQUIRE(result.rejected_steps == 0U);
    // One tetrahedron has four degree-two volume quadrature points. The
    // spatial callbacks are evaluated once per study, not once per trial.
    NOVA_REQUIRE(conductivity_samples == 4U);
    NOVA_REQUIRE(capacity_samples == 4U);
    NOVA_REQUIRE(source_samples == 4U);
}
