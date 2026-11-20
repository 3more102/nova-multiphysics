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


NOVA_TEST(transient_heat_3d_time_dependent_uniform_source_uses_scheme_quadrature) {
    auto p = unit_heat();
    p.diffusion.source = {};
    p.diffusion.cell_source.clear();
    p.cell_volumetric_capacity = {2.0};
    p.initial_temperature = {300.0, 300.0, 300.0, 300.0};
    std::size_t samples = 0U;
    p.time_dependent_source = [&](double, double, double, double time) {
        ++samples;
        return time;
    };

    TransientHeat3DOptions options;
    options.start_time = 2.0;
    options.time_step = 0.1;
    options.steps = 10U;
    options.solver_options.relative_tolerance = 1e-13;
    options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;

    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                              TransientHeat3DScheme::CrankNicolson}) {
        samples = 0U;
        options.scheme = scheme;
        const auto result = solve_transient_heat_3d(p, options);
        const bool backward_euler = scheme == TransientHeat3DScheme::BackwardEuler;
        const double expected_rise = backward_euler ? 1.275 : 1.25;
        const double expected_source_energy = backward_euler ? 17.0 / 40.0 : 5.0 / 12.0;
        NOVA_REQUIRE_NEAR(result.final_time, 3.0, 1e-14);
        for (const double temperature : result.temperature) {
            NOVA_REQUIRE_NEAR(temperature, 300.0 + expected_rise, 2e-10);
        }
        double total_source_energy = 0.0;
        for (const auto& diagnostic : result.diagnostics) {
            total_source_energy += diagnostic.source_energy;
            NOVA_REQUIRE_NEAR(diagnostic.balance_residual, 0.0, 1e-10);
        }
        NOVA_REQUIRE_NEAR(total_source_energy, expected_source_energy, 1e-12);
        NOVA_REQUIRE_NEAR(result.diagnostics.back().heat_content,
                          100.0 + expected_source_energy, 2e-10);
        NOVA_REQUIRE(samples == (backward_euler ? 40U : 44U));
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
        bad.time_dependent_source = [](double, double, double, double) {
            return std::numeric_limits<double>::quiet_NaN();
        };
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


NOVA_TEST(transient_heat_3d_cellwise_tensor_capacity_source_patch) {
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 1, 1, 1);
    auto p = unit_heat();
    p.diffusion.mesh = box.mesh;
    p.diffusion.cell_diffusivity.clear();
    p.diffusion.source = {};
    p.cell_volumetric_capacity.clear();
    p.diffusion.cell_source.clear();
    for (std::size_t i = 0; i < box.mesh.elements.size(); ++i) {
        const double capacity = 2.0 + static_cast<double>(i);
        p.cell_volumetric_capacity.push_back(capacity);
        p.diffusion.cell_source.push_back(4.0 * capacity);
        p.diffusion.cell_diffusivity_tensors.push_back({
            1.0 + static_cast<double>(i), 2.0 + static_cast<double>(i),
            3.0 + static_cast<double>(i), 0.1, 0.0, 0.05
        });
    }
    p.initial_temperature.assign(box.mesh.nodes.size(), 300.0);
    TransientHeat3DOptions options;
    options.time_step = 0.25;
    options.solver_options.relative_tolerance = 1e-12;
    options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                               TransientHeat3DScheme::CrankNicolson}) {
        options.scheme = scheme;
        const auto solved = solve_transient_heat_3d(p, options);
        for (const double temperature : solved.temperature) {
            NOVA_REQUIRE_NEAR(temperature, 301.0, 2e-9);
        }
        NOVA_REQUIRE_NEAR(solved.diagnostics[0].source_energy, 4.5, 2e-11);
        NOVA_REQUIRE_NEAR(solved.diagnostics[0].heat_content, 1354.5, 2e-8);
        NOVA_REQUIRE_NEAR(solved.diagnostics[0].balance_residual, 0.0, 1e-8);
    }
}

NOVA_TEST(transient_heat_3d_missing_spatial_and_cell_sources_means_zero_generation) {
    auto p = unit_heat();
    p.diffusion.source = {};
    p.diffusion.cell_source.clear();
    p.initial_temperature = {291.0, 291.0, 291.0, 291.0};
    TransientHeat3DOptions options;
    options.steps = 3;
    for (auto scheme : {TransientHeat3DScheme::BackwardEuler,
                        TransientHeat3DScheme::CrankNicolson}) {
        options.scheme = scheme;
        const auto result = solve_transient_heat_3d(p, options);
        for (double value : result.temperature) {
            NOVA_REQUIRE_NEAR(value, 291.0, 1e-11);
        }
        for (const auto& budget : result.diagnostics) {
            NOVA_REQUIRE_NEAR(budget.source_energy, 0.0, 1e-15);
            NOVA_REQUIRE_NEAR(budget.balance_residual, 0.0, 1e-10);
        }
    }
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


NOVA_TEST(transient_heat_3d_time_dependent_dirichlet_node_tracks_physical_time) {
    auto p = unit_heat();
    p.diffusion.dirichlet_constraints = {{0U, 0.0}, {1U, 0.0}, {2U, 0.0}};
    p.time_dependent_dirichlet_constraints.push_back({
        3U, [](double, double, double, double time) { return 2.0 * time; }
    });

    TransientHeat3DOptions options;
    options.start_time = 2.0;
    options.time_step = 0.1;
    options.steps = 5U;
    options.solver_options.relative_tolerance = 1e-13;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::Jacobi;

    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                              TransientHeat3DScheme::CrankNicolson}) {
        options.scheme = scheme;
        const auto result = solve_transient_heat_3d(p, options);
        NOVA_REQUIRE_NEAR(result.final_time, 2.5, 1e-14);
        for (std::size_t node = 0; node < 3U; ++node) {
            NOVA_REQUIRE_NEAR(result.temperature[node], 0.0, 1e-14);
        }
        NOVA_REQUIRE_NEAR(result.temperature[3], 5.0, 1e-12);
        NOVA_REQUIRE_NEAR(result.diagnostics.back().heat_content,
                          5.0 / 24.0, 1e-11);
        double imposed_energy = 0.0;
        for (const auto& diagnostic : result.diagnostics) {
            imposed_energy += diagnostic.dirichlet_energy;
            NOVA_REQUIRE_NEAR(diagnostic.balance_residual, 0.0, 2e-10);
        }
        NOVA_REQUIRE_NEAR(imposed_energy, 1.0 / 24.0, 2e-10);
    }
}

NOVA_TEST(transient_heat_3d_time_dependent_neumann_uses_scheme_time_quadrature) {
    auto p = unit_heat();
    p.time_dependent_neumann_faces.push_back({
        {0U, 1U, 2U},
        [](double, double, double, double time) { return time; }
    });

    TransientHeat3DOptions options;
    options.start_time = 0.0;
    options.time_step = 0.1;
    options.steps = 10U;
    options.solver_options.relative_tolerance = 1e-13;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::Jacobi;

    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                              TransientHeat3DScheme::CrankNicolson}) {
        options.scheme = scheme;
        const auto result = solve_transient_heat_3d(p, options);
        const double expected_energy =
            scheme == TransientHeat3DScheme::BackwardEuler ? -0.275 : -0.25;
        double total_neumann_energy = 0.0;
        for (const auto& diagnostic : result.diagnostics) {
            total_neumann_energy += diagnostic.neumann_energy;
            NOVA_REQUIRE_NEAR(diagnostic.balance_residual, 0.0, 2e-10);
        }
        NOVA_REQUIRE_NEAR(total_neumann_energy, expected_energy, 2e-12);
        NOVA_REQUIRE_NEAR(result.diagnostics.back().heat_content,
                          expected_energy, 2e-10);
    }
}

NOVA_TEST(transient_heat_3d_rejects_invalid_time_dependent_boundary_data) {
    auto bad_dirichlet = unit_heat();
    bad_dirichlet.time_dependent_dirichlet_constraints.push_back({
        3U, [](double, double, double, double) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    });
    NOVA_REQUIRE(rejects([&] {
        (void)solve_transient_heat_3d(bad_dirichlet);
    }));

    auto bad_neumann = unit_heat();
    bad_neumann.time_dependent_neumann_faces.push_back({
        {0U, 1U, 2U},
        [](double, double, double, double) {
            return std::numeric_limits<double>::infinity();
        }
    });
    NOVA_REQUIRE(rejects([&] {
        (void)solve_transient_heat_3d(bad_neumann);
    }));

    auto missing_callback = unit_heat();
    missing_callback.time_dependent_dirichlet_faces.push_back({{0U, 1U, 2U}, {}});
    NOVA_REQUIRE(rejects([&] {
        (void)solve_transient_heat_3d(missing_callback);
    }));
}
