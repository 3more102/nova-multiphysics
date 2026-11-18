#include "test_framework.hpp"

#include "nova/fem/adaptive_heat_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

nova::fem::TransientHeat3DProblem unit_tetrahedron() {
    using namespace nova::fem;
    TransientHeat3DProblem problem;
    problem.mesh.nodes = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
    problem.mesh.elements = {{0,1,2,3}};
    problem.volumetric_capacity = [](double, double, double) { return 1.0; };
    problem.conductivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double, double) { return 0.0; };
    problem.initial_temperature = [](double, double, double) { return 0.0; };
    return problem;
}

} // namespace

NOVA_TEST(adaptive_heat_3d_preserves_disconnected_component_heat_and_reuses_material_operators) {
    using namespace nova::fem;
    auto problem = unit_tetrahedron();
    problem.mesh.nodes.insert(problem.mesh.nodes.end(), {
        {2,0,0}, {3,0,0}, {2,1,0}, {2,0,1}
    });
    problem.mesh.elements.push_back({4,5,6,7});
    int capacity_samples = 0;
    int conductivity_samples = 0;
    problem.volumetric_capacity = [&](double, double, double) {
        ++capacity_samples;
        return 1.0;
    };
    problem.conductivity = [&](double, double, double) {
        ++conductivity_samples;
        return 1.0;
    };
    problem.initial_temperature = [](double x, double, double) {
        return x < 2.0 ? 3.0 : 5.0;
    };

    AdaptiveHeat3DOptions options;
    options.end_time = 0.5;
    options.initial_step = 0.2;
    options.max_step = 0.2;
    options.absolute_tolerance = 1.0e-9;
    options.relative_tolerance = 0.0;
    const auto first = solve_adaptive_heat_3d(problem, options);

    NOVA_REQUIRE(first.accepted_steps >= 2U);
    NOVA_REQUIRE(first.rejected_steps == 0U);
    NOVA_REQUIRE(first.diagnostics.size() == first.accepted_steps);
    NOVA_REQUIRE(first.substep_solvers.size() == 2U * first.accepted_steps);
    for (std::size_t i = 0; i < first.temperature.size(); ++i) {
        NOVA_REQUIRE_NEAR(first.temperature[i], i < 4U ? 3.0 : 5.0, 1.0e-11);
    }
    for (const auto& diagnostic : first.diagnostics) {
        NOVA_REQUIRE_NEAR(diagnostic.heat_content, 8.0 / 6.0, 1.0e-11);
        NOVA_REQUIRE_NEAR(diagnostic.balance_residual, 0.0, 1.0e-11);
    }

    // An exact DOF restart must preserve each insulated component independently.
    problem.initial_dof_values = first.temperature;
    problem.initial_temperature = {};
    options.start_time = first.final_time;
    options.end_time = 1.0;
    const auto second = solve_adaptive_heat_3d(problem, options);
    NOVA_REQUIRE_NEAR(second.final_time, 1.0, 1.0e-14);
    for (std::size_t i = 0; i < second.temperature.size(); ++i) {
        NOVA_REQUIRE_NEAR(second.temperature[i], i < 4U ? 3.0 : 5.0, 1.0e-11);
    }
    NOVA_REQUIRE(capacity_samples == 16);
    NOVA_REQUIRE(conductivity_samples == 16);
}

NOVA_TEST(adaptive_heat_3d_backward_euler_rejects_large_steps_for_time_dependent_source) {
    using namespace nova::fem;
    auto problem = unit_tetrahedron();
    problem.source = [](double, double, double, double time) { return 2.0 * time; };

    int capacity_samples = 0;
    int conductivity_samples = 0;
    problem.volumetric_capacity = [&](double, double, double) {
        ++capacity_samples;
        return 1.0;
    };
    problem.conductivity = [&](double, double, double) {
        ++conductivity_samples;
        return 1.0;
    };

    AdaptiveHeat3DOptions options;
    options.scheme = TransientTimeScheme::BackwardEuler;
    options.end_time = 0.8;
    options.initial_step = 0.4;
    options.max_step = 0.4;
    options.min_step = 1.0e-5;
    options.absolute_tolerance = 2.0e-4;
    options.relative_tolerance = 0.0;
    const auto result = solve_adaptive_heat_3d(problem, options);

    NOVA_REQUIRE(result.accepted_steps > 2U);
    NOVA_REQUIRE(result.rejected_steps > 0U);
    NOVA_REQUIRE(result.attempts.size() == result.accepted_steps + result.rejected_steps);
    NOVA_REQUIRE_NEAR(result.final_time, 0.8, 1.0e-14);
    for (const double temperature : result.temperature) {
        NOVA_REQUIRE_NEAR(temperature, 0.64, 0.03);
    }
    for (const auto& diagnostic : result.diagnostics) {
        NOVA_REQUIRE_NEAR(diagnostic.balance_residual, 0.0, 2.0e-9);
    }
    // Four tetrahedron quadrature points per material callback, once per study.
    NOVA_REQUIRE(capacity_samples == 4);
    NOVA_REQUIRE(conductivity_samples == 4);
}

NOVA_TEST(adaptive_heat_3d_crank_nicolson_preserves_nonzero_time_restart_and_face_data) {
    using namespace nova::fem;
    auto problem = unit_tetrahedron();
    constexpr double start = 0.4;
    problem.initial_dof_values = {start, 1.0 + start, start, start};
    problem.initial_temperature = {};
    problem.source = [](double, double, double, double) { return 1.0; };
    problem.dirichlet_faces.push_back({
        {0,2,3},
        [](double, double, double, double time) { return time; }
    });
    problem.neumann_faces.push_back({
        {1,2,3},
        [](double, double, double, double) { return -1.0 / std::sqrt(3.0); }
    });

    AdaptiveHeat3DOptions options;
    options.scheme = TransientTimeScheme::CrankNicolson;
    options.start_time = start;
    options.end_time = 0.9;
    options.initial_step = 0.25;
    options.max_step = 0.25;
    options.absolute_tolerance = 1.0e-8;
    options.relative_tolerance = 0.0;
    const auto result = solve_adaptive_heat_3d(problem, options);

    NOVA_REQUIRE_NEAR(result.final_time, 0.9, 1.0e-14);
    for (std::size_t i = 0; i < result.temperature.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.temperature[i], problem.mesh.nodes[i].x + 0.9, 2.0e-10);
    }
    NOVA_REQUIRE(result.diagnostics.size() == result.accepted_steps);
    for (const auto& diagnostic : result.diagnostics) {
        NOVA_REQUIRE_NEAR(diagnostic.balance_residual, 0.0, 2.0e-10);
    }
}

NOVA_TEST(adaptive_heat_3d_crank_nicolson_tightening_improves_against_refined_fixed_step) {
    using namespace nova::fem;
    auto problem = unit_tetrahedron();
    problem.initial_temperature = [](double x, double, double) { return x; };
    problem.dirichlet_faces.push_back({
        {0,2,3},
        [](double, double, double, double) { return 0.0; }
    });

    TransientHeat3DOptions reference_options;
    reference_options.scheme = TransientTimeScheme::CrankNicolson;
    reference_options.time_step = 0.0005;
    reference_options.steps = 400U;
    reference_options.solver_options.relative_tolerance = 1.0e-13;
    const auto reference = solve_transient_heat_3d(problem, reference_options);

    AdaptiveHeat3DOptions loose_options;
    loose_options.scheme = TransientTimeScheme::CrankNicolson;
    loose_options.end_time = 0.2;
    loose_options.initial_step = 0.2;
    loose_options.max_step = 0.2;
    loose_options.min_step = 1.0e-6;
    loose_options.absolute_tolerance = 1.0e-3;
    loose_options.relative_tolerance = 0.0;
    const auto loose = solve_adaptive_heat_3d(problem, loose_options);

    auto strict_options = loose_options;
    strict_options.absolute_tolerance = 1.0e-6;
    const auto strict = solve_adaptive_heat_3d(problem, strict_options);

    const double loose_error = std::abs(loose.temperature[1] - reference.temperature[1]);
    const double strict_error = std::abs(strict.temperature[1] - reference.temperature[1]);
    NOVA_REQUIRE(strict.accepted_steps > loose.accepted_steps);
    NOVA_REQUIRE(strict.rejected_steps > 0U);
    NOVA_REQUIRE(strict_error < loose_error);
    NOVA_REQUIRE_NEAR(strict.final_time, 0.2, 1.0e-14);
}

NOVA_TEST(adaptive_heat_3d_rejects_invalid_controls_restart_and_attempt_exhaustion) {
    using namespace nova::fem;
    auto problem = unit_tetrahedron();
    AdaptiveHeat3DOptions options;
    options.end_time = options.start_time;
    bool invalid_time = false;
    try {
        (void)solve_adaptive_heat_3d(problem, options);
    } catch (const std::invalid_argument&) {
        invalid_time = true;
    }
    NOVA_REQUIRE(invalid_time);

    options.end_time = 1.0;
    options.absolute_tolerance = 0.0;
    bool invalid_tolerance = false;
    try {
        (void)solve_adaptive_heat_3d(problem, options);
    } catch (const std::invalid_argument&) {
        invalid_tolerance = true;
    }
    NOVA_REQUIRE(invalid_tolerance);

    options.absolute_tolerance = 1.0e-14;
    options.scheme = TransientTimeScheme::BackwardEuler;
    options.initial_step = 0.8;
    options.max_step = 0.8;
    options.max_attempts = 1U;
    problem.source = [](double, double, double, double time) { return 2.0 * time; };
    bool attempts_exhausted = false;
    try {
        (void)solve_adaptive_heat_3d(problem, options);
    } catch (const std::runtime_error&) {
        attempts_exhausted = true;
    }
    NOVA_REQUIRE(attempts_exhausted);

    problem = unit_tetrahedron();
    problem.initial_dof_values = {0.0, 1.0, 0.0, std::numeric_limits<double>::infinity()};
    bool invalid_restart = false;
    try {
        (void)solve_adaptive_heat_3d(problem);
    } catch (const std::invalid_argument&) {
        invalid_restart = true;
    }
    NOVA_REQUIRE(invalid_restart);
}
