#include "test_framework.hpp"

#include "nova/fem/adaptive_heat_3d.hpp"
#include "nova/fem/diffusion_3d.hpp"
#include "nova/linalg/bicgstab.hpp"
#include "nova/linalg/conjugate_gradient.hpp"
#include "nova/linalg/gmres.hpp"
#include "nova/linalg/incomplete_cholesky.hpp"

#include <array>
#include <limits>

namespace {
using namespace nova::linalg;

template<class F>
bool rejects(F&& action) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

SolverOptions ic_options() {
    SolverOptions options;
    options.preconditioner = PreconditionerType::IncompleteCholesky0;
    options.relative_tolerance = 1e-12;
    return options;
}
} // namespace

NOVA_TEST(ic0_gives_the_known_incomplete_inverse_on_a_four_cycle) {
    const auto matrix = CsrMatrix::from_triplets(4, 4, {
        {0,0,4}, {1,1,4}, {2,2,4}, {3,3,4}, {0,1,-1}, {1,0,-1},
        {1,2,-1}, {2,1,-1}, {2,3,-1}, {3,2,-1}, {0,3,-1}, {3,0,-1}});
    const IncompleteCholesky0 factor(matrix);
    NOVA_REQUIRE(factor.dimension() == 4 && factor.stored_entries() == 8);
    // LL^T matches A on its original pattern. L10 = L30 = -1/2,
    // giving (LL^T)13 = (LL^T)31 = 1/4, whereas A13 = A31 = 0.
    for (const auto& exact : std::array{std::array{1.,2.,-3.,4.}, std::array{0.,-4.,2.,1.}}) {
        auto rhs = matrix.multiply(exact);
        rhs[1] += exact[3]/4;
        rhs[3] += exact[1]/4;
        const auto result = factor.apply(rhs);
        for (std::size_t i = 0; i < 4; ++i) NOVA_REQUIRE_NEAR(result[i], exact[i], 3e-15);
    }
}

NOVA_TEST(ic0_solves_a_full_spd_pattern_in_one_cg_iteration_and_reports_true_residual) {
    // A = L L^T for an integer lower triangular L with nonzero row intersections.
    const auto matrix = CsrMatrix::from_triplets(4, 4, {
        {0,0,4}, {0,1,2}, {0,2,-2}, {0,3,4}, {1,0,2}, {1,1,5}, {1,2,1}, {1,3,4},
        {2,0,-2}, {2,1,1}, {2,2,11}, {2,3,2}, {3,0,4}, {3,1,4}, {3,2,2}, {3,3,10}});
    const std::array exact{1.,-2.,3.,0.25};
    const auto rhs = matrix.multiply(exact);
    const auto result = solve_conjugate_gradient(matrix, rhs, ic_options(), std::array{0.,1.,0.,0.});
    NOVA_REQUIRE(result.converged && result.iterations == 1);
    NOVA_REQUIRE(result.message.find("IC(0)") != std::string::npos);
    const auto product = matrix.multiply(result.solution);
    double true_norm = 0;
    for (std::size_t i = 0; i < exact.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.solution[i], exact[i], 4e-14);
        true_norm = std::hypot(true_norm, rhs[i]-product[i]);
    }
    NOVA_REQUIRE_NEAR(result.residual_norm, true_norm, 1e-26);
}

NOVA_TEST(ic0_accepts_well_conditioned_correlations_across_material_units) {
    for (double s : {1e-240, 1.0, 1e240}) {
        const auto matrix = CsrMatrix::from_triplets(2,2,{{0,0,4*s},{0,1,s},{1,0,s},{1,1,3*s}});
        const auto result = IncompleteCholesky0(matrix).apply(std::array{2*s,-5*s});
        NOVA_REQUIRE_NEAR(result[0], 1, 2e-15);
        NOVA_REQUIRE_NEAR(result[1], -2, 2e-15);
    }
    const auto diagonal = CsrMatrix::from_triplets(2,2,{{0,0,1e-240},{1,1,1e240}});
    const auto scaled = IncompleteCholesky0(diagonal).apply(std::array{1e-240,-2e240});
    NOVA_REQUIRE_NEAR(scaled[0], 1, 2e-15);
    NOVA_REQUIRE_NEAR(scaled[1], -2, 2e-15);
}

NOVA_TEST(ic0_rejects_asymmetry_unresolved_pivots_and_invalid_applications) {
    NOVA_REQUIRE(rejects([] { (void)IncompleteCholesky0(CsrMatrix::from_triplets(2,3,{})); }));
    NOVA_REQUIRE(rejects([] { (void)IncompleteCholesky0(CsrMatrix::from_triplets(2,2,{{0,0,1}})); }));
    for (const auto& entries : std::vector<std::vector<Triplet>>{
             {{0,0,2},{1,1,2},{0,1,0.5}}, {{0,0,2},{1,1,2},{1,0,0.5}},
             {{0,0,2},{1,1,2},{0,1,0.5},{1,0,0.6}},
             {{0,0,2},{1,1,2},{0,1,2},{1,0,2}}}) {
        NOVA_REQUIRE(rejects([&] { (void)IncompleteCholesky0(CsrMatrix::from_triplets(2,2,entries)); }));
    }
    const auto identity = CsrMatrix::from_triplets(1,1,{{0,0,1}});
    const IncompleteCholesky0 factor(identity);
    NOVA_REQUIRE(rejects([&] { (void)factor.apply({}); }));
    NOVA_REQUIRE(rejects([&] { (void)factor.apply(std::array{std::numeric_limits<double>::infinity()}); }));
    GmresOptions gmres;
    gmres.solver = ic_options();
    NOVA_REQUIRE(rejects([&] { (void)solve_gmres(identity, std::array{1.}, gmres); }));
    NOVA_REQUIRE(rejects([&] { (void)solve_bicgstab(identity, std::array{1.}, ic_options()); }));
}

NOVA_TEST(ic0_reports_breakdown_on_an_spd_matrix_whose_zero_fill_pivot_is_negative) {
    // Leading principal minors 4, 27, 99 and 54 are positive. Omission of
    // the fill at (3,1) nevertheless gives a negative final IC(0) pivot.
    const auto matrix = CsrMatrix::from_triplets(4,4,{
        {0,0,4}, {1,1,9}, {2,2,9}, {3,3,9}, {0,1,3}, {1,0,3},
        {1,2,-6}, {2,1,-6}, {2,3,-4}, {3,2,-4}, {0,3,-5}, {3,0,-5}});
    NOVA_REQUIRE(rejects([&] { (void)solve_conjugate_gradient(matrix, std::array{1.,1.,1.,1.}, ic_options()); }));
}

NOVA_TEST(ic0_accelerates_a_tetrahedral_diffusion_solution_at_the_same_tolerance) {
    using namespace nova::fem;
    auto box = make_box_mesh_3d(0,1,0,1,0,1, 6,6,6);
    ScalarDiffusion3DProblem p;
    p.mesh = box.mesh;
    p.diffusivity_tensor = [](double,double,double) { return SymmetricTensor3D{2,3,4,0.25,-0.2,0.15}; };
    p.source = [](double,double,double) { return 0.0; };
    const auto exact = [](double x,double y,double z) { return 1+x-2*y+z; };
    for (const auto& face : box.boundary_faces) p.dirichlet_faces.push_back({face, exact});
    auto options = ic_options();
    const auto ic = solve_scalar_diffusion_3d(p, options);
    options.preconditioner = PreconditionerType::Jacobi;
    const auto jacobi = solve_scalar_diffusion_3d(p, options);
    NOVA_REQUIRE(ic.linear_solver.converged && jacobi.linear_solver.converged);
    NOVA_REQUIRE(ic.linear_solver.iterations < jacobi.linear_solver.iterations);
    NOVA_REQUIRE(ic.linear_solver.residual_norm <= options.relative_tolerance * ic.linear_solver.initial_residual_norm);
    for (std::size_t i = 0; i < p.mesh.nodes.size(); ++i) {
        const auto& n = p.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(ic.values[i], exact(n.x,n.y,n.z), 1e-10);
        NOVA_REQUIRE_NEAR(ic.values[i], jacobi.values[i], 1e-10);
    }
}

NOVA_TEST(ic0_and_jacobi_agree_on_a_transient_tensor_heat_history) {
    using namespace nova::fem;
    const auto box = make_box_mesh_3d(0,1,0,1,0,1, 3,3,3);
    TransientHeat3DProblem p;
    p.mesh = box.mesh;
    p.volumetric_capacity = [](double,double,double) { return 1.0; };
    p.conductivity_tensor = [](double x,double,double) { return SymmetricTensor3D{2+x,3,1,0.3,0.1,-0.2}; };
    p.initial_temperature = [](double x,double y,double z) { return 64*x*(1-x)*y*(1-y)*z*(1-z); };
    p.source = [](double,double,double,double) { return 0.0; };
    for (const auto& face : box.boundary_faces) p.dirichlet_faces.push_back({face, p.source});
    AdaptiveHeat3DOptions options;
    options.end_time = options.max_step = options.initial_step = 0.025;
    options.absolute_tolerance = 1e-5;
    options.relative_tolerance = 0;
    options.solver_options = ic_options();
    const auto ic = solve_adaptive_heat_3d(p, options);
    options.solver_options.preconditioner = PreconditionerType::Jacobi;
    const auto jacobi = solve_adaptive_heat_3d(p, options);
    NOVA_REQUIRE(ic.rejected_steps > 0 && ic.accepted_steps > 1);
    NOVA_REQUIRE(ic.accepted_steps == jacobi.accepted_steps && ic.rejected_steps == jacobi.rejected_steps);
    for (std::size_t i = 0; i < p.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(ic.temperature[i], jacobi.temperature[i], 2e-10);
    }
    for (const auto& d : ic.diagnostics) NOVA_REQUIRE_NEAR(d.balance_residual, 0, 2e-11);
    for (const auto& solve : ic.substep_solvers) {
        NOVA_REQUIRE(solve.converged && solve.message.find("IC(0)") != std::string::npos);
    }
}

NOVA_TEST(ic0_cache_reuses_an_exact_matrix_and_rebuilds_after_changes) {
    const auto matrix = CsrMatrix::from_triplets(3, 3, {
        {0,0,4}, {1,1,5}, {2,2,6},
        {0,1,-1}, {1,0,-1}, {1,2,-1}, {2,1,-1}});
    const auto changed = CsrMatrix::from_triplets(3, 3, {
        {0,0,4}, {1,1,5.5}, {2,2,6},
        {0,1,-1}, {1,0,-1}, {1,2,-1}, {2,1,-1}});
    const auto diagonal = CsrMatrix::from_triplets(3, 3, {
        {0,0,4}, {1,1,5}, {2,2,6}});
    const std::array solution_one{1.0, -2.0, 3.0};
    const std::array solution_two{-4.0, 3.0, 1.0};
    const auto opts = ic_options();
    CachedConjugateGradient cache;

    const auto first = cache.solve(matrix, matrix.multiply(solution_one), opts);
    NOVA_REQUIRE(first.converged && cache.factor_builds() == 1);
    const auto second = cache.solve(matrix, matrix.multiply(solution_two), opts);
    NOVA_REQUIRE(second.converged && cache.factor_builds() == 1);
    for (std::size_t i = 0; i < 3; ++i) {
        NOVA_REQUIRE_NEAR(first.solution[i], solution_one[i], 1e-11);
        NOVA_REQUIRE_NEAR(second.solution[i], solution_two[i], 1e-11);
    }

    const auto third = cache.solve(changed, changed.multiply(solution_one), opts);
    NOVA_REQUIRE(third.converged && cache.factor_builds() == 2);
    const auto fourth = cache.solve(diagonal, diagonal.multiply(solution_two), opts);
    NOVA_REQUIRE(fourth.converged && cache.factor_builds() == 3);
    for (std::size_t i = 0; i < 3; ++i) {
        NOVA_REQUIRE_NEAR(third.solution[i], solution_one[i], 1e-11);
        NOVA_REQUIRE_NEAR(fourth.solution[i], solution_two[i], 1e-11);
    }
    // A different solver choice cannot accidentally use the stored IC(0) factor.
    auto jacobi = opts;
    jacobi.preconditioner = PreconditionerType::Jacobi;
    NOVA_REQUIRE(cache.solve(matrix, matrix.multiply(solution_one), jacobi).converged);
    NOVA_REQUIRE(cache.factor_builds() == 3);
    NOVA_REQUIRE(cache.solve(matrix, matrix.multiply(solution_one), opts).converged);
    NOVA_REQUIRE(cache.factor_builds() == 4);
}

NOVA_TEST(ic0_fixed_step_heat_agrees_with_jacobi_for_changing_robin_matrix) {
    using namespace nova::fem;
    auto p = TransientHeat3DProblem{};
    p.mesh = make_box_mesh_3d(0,1,0,1,0,1, 2,2,2).mesh;
    p.volumetric_capacity = [](double,double,double) { return 1.0; };
    p.conductivity = [](double,double,double) { return 2.0; };
    p.initial_temperature = [](double x,double y,double) { return x+y; };
    p.source = [](double,double,double,double) { return 0.0; };
    const auto box = make_box_mesh_3d(0,1,0,1,0,1, 2,2,2);
    for (const auto& face : box.boundary_faces) {
        p.robin_faces.push_back({face,
            [](double,double,double,double time) { return 1.0 + time; },
            [](double,double,double,double time) { return 3.0 - time; }});
    }
    TransientHeat3DOptions o;
    o.steps = 4;
    o.time_step = 0.125;
    o.solver_options = ic_options();
    const auto ic = solve_transient_heat_3d(p, o);
    o.solver_options.preconditioner = PreconditionerType::Jacobi;
    const auto jacobi = solve_transient_heat_3d(p, o);
    for (std::size_t i = 0; i < ic.temperature.size(); ++i) {
        NOVA_REQUIRE_NEAR(ic.temperature[i], jacobi.temperature[i], 1e-10);
    }
    for (const auto& solver : ic.linear_solvers) {
        NOVA_REQUIRE(solver.converged);
    }
}


namespace {
nova::fem::TransientHeat3DProblem uniform_insulated_heat_tet() {
    nova::fem::TransientHeat3DProblem problem;
    problem.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    problem.volumetric_capacity = [](double, double, double) { return 1.0; };
    problem.conductivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double, double) { return 2.0; };
    problem.initial_temperature = [](double, double, double) { return 0.0; };
    return problem;
}
} // namespace

NOVA_TEST(adaptive_3d_ic0_reuses_full_and_half_factors_across_macro_steps) {
    using namespace nova::fem;
    auto problem = uniform_insulated_heat_tet();
    AdaptiveHeat3DOptions options;
    options.end_time = 0.5;
    options.initial_step = options.min_step = options.max_step = 0.125;
    options.absolute_tolerance = 1e-6;
    options.relative_tolerance = 0.0;
    options.solver_options = ic_options();
    for (auto scheme : {TransientTimeScheme::BackwardEuler, TransientTimeScheme::CrankNicolson}) {
        options.scheme = scheme;
        const auto result = solve_adaptive_heat_3d(problem, options);
        NOVA_REQUIRE(result.accepted_steps == 4 && result.rejected_steps == 0);
        NOVA_REQUIRE(result.substep_solvers.size() == 8);
        NOVA_REQUIRE(result.full_step_ic0_factor_builds == 1);
        NOVA_REQUIRE(result.half_step_ic0_factor_builds == 1);
        for (const auto value : result.temperature) {
            NOVA_REQUIRE_NEAR(value, 1.0, 2e-11);
        }
        for (const auto& solver : result.substep_solvers) {
            NOVA_REQUIRE(solver.converged);
        }
        auto jacobi = options;
        jacobi.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
        const auto baseline = solve_adaptive_heat_3d(problem, jacobi);
        NOVA_REQUIRE(baseline.full_step_ic0_factor_builds == 0);
        NOVA_REQUIRE(baseline.half_step_ic0_factor_builds == 0);
        for (std::size_t i = 0; i < result.temperature.size(); ++i) {
            NOVA_REQUIRE_NEAR(result.temperature[i], baseline.temperature[i], 2e-11);
        }
    }
}

NOVA_TEST(adaptive_3d_ic0_rebuilds_for_different_full_and_half_time_steps) {
    using namespace nova::fem;
    const auto problem = uniform_insulated_heat_tet();
    AdaptiveHeat3DOptions options;
    options.end_time = 0.375;
    options.initial_step = options.max_step = 0.25;
    options.min_step = 0.0625;
    options.absolute_tolerance = 1e-6;
    options.relative_tolerance = 0.0;
    options.scheme = TransientTimeScheme::BackwardEuler;
    options.solver_options = ic_options();
    const auto result = solve_adaptive_heat_3d(problem, options);
    NOVA_REQUIRE(result.accepted_steps == 2 && result.rejected_steps == 0);
    NOVA_REQUIRE(result.full_step_ic0_factor_builds == 2);
    NOVA_REQUIRE(result.half_step_ic0_factor_builds == 2);
    for (const auto value : result.temperature) {
        NOVA_REQUIRE_NEAR(value, 0.75, 2e-11);
    }
}

NOVA_TEST(adaptive_3d_ic0_rebuilds_when_robin_transfer_changes_with_time) {
    using namespace nova::fem;
    auto problem = uniform_insulated_heat_tet();
    problem.robin_faces.push_back({{1,2,3},
        [](double, double, double, double time) { return 1.0 + time; },
        [](double, double, double, double) { return 0.0; }});
    AdaptiveHeat3DOptions options;
    options.end_time = 0.25;
    options.initial_step = options.min_step = options.max_step = 0.125;
    options.absolute_tolerance = 1.0;
    options.relative_tolerance = 0.0;
    options.scheme = TransientTimeScheme::BackwardEuler;
    options.solver_options = ic_options();
    const auto ic = solve_adaptive_heat_3d(problem, options);
    NOVA_REQUIRE(ic.accepted_steps == 2 && ic.rejected_steps == 0);
    NOVA_REQUIRE(ic.full_step_ic0_factor_builds == 2);
    NOVA_REQUIRE(ic.half_step_ic0_factor_builds == 4);
    auto jacobi = options;
    jacobi.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    const auto baseline = solve_adaptive_heat_3d(problem, jacobi);
    for (std::size_t i = 0; i < ic.temperature.size(); ++i) {
        NOVA_REQUIRE_NEAR(ic.temperature[i], baseline.temperature[i], 2e-10);
    }
    for (const auto& diagnostic : ic.diagnostics) {
        NOVA_REQUIRE_NEAR(diagnostic.balance_residual, 0.0, 2e-11);
    }
}
