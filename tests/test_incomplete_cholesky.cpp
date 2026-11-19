#include "test_framework.hpp"

#include "nova/fem/transient_heat_3d.hpp"
#include "nova/fem/diffusion_3d.hpp"
#include "nova/linalg/bicgstab.hpp"
#include "nova/linalg/gmres.hpp"
#include "nova/linalg/incomplete_cholesky.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

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

NOVA_TEST(ic0_and_jacobi_agree_on_fixed_step_transient_heat_histories) {
    using namespace nova::fem;
    const auto box = make_box_mesh_3d(0,1,0,1,0,1, 3,3,3);

    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                              TransientHeat3DScheme::CrankNicolson}) {
        TransientHeat3DProblem p;
        p.diffusion.mesh = box.mesh;
        p.diffusion.diffusivity_tensor = [](double x,double,double) {
            return SymmetricTensor3D{2+x,3,1,0.3,0.1,-0.2};
        };
        p.volumetric_capacity = [](double,double,double) { return 1.0; };
        p.diffusion.source = [](double,double,double) { return 0.0; };
        p.initial_temperature.resize(p.diffusion.mesh.nodes.size());
        for (std::size_t i = 0; i < p.diffusion.mesh.nodes.size(); ++i) {
            const auto& node = p.diffusion.mesh.nodes[i];
            p.initial_temperature[i] =
                64*node.x*(1-node.x)*node.y*(1-node.y)*node.z*(1-node.z);
        }
        for (const auto& face : box.boundary_faces) {
            p.diffusion.dirichlet_faces.push_back({face,
                [](double,double,double) { return 0.0; }});
        }

        TransientHeat3DOptions options;
        options.scheme = scheme;
        options.start_time = 0.5;
        options.time_step = 0.025;
        options.steps = 4;
        options.solver_options = ic_options();
        const auto ic = solve_transient_heat_3d(p, options);
        options.solver_options.preconditioner = PreconditionerType::Jacobi;
        const auto jacobi = solve_transient_heat_3d(p, options);

        NOVA_REQUIRE(ic.linear_solvers.size() == options.steps);
        NOVA_REQUIRE(ic.diagnostics.size() == options.steps);
        NOVA_REQUIRE(ic.final_time == jacobi.final_time);
        NOVA_REQUIRE_NEAR(ic.final_time, 0.6, 1e-12);
        for (std::size_t i = 0; i < ic.temperature.size(); ++i) {
            NOVA_REQUIRE_NEAR(ic.temperature[i], jacobi.temperature[i], 2e-10);
        }
        for (std::size_t step = 0; step < options.steps; ++step) {
            const auto& ic_solve = ic.linear_solvers[step];
            NOVA_REQUIRE(ic_solve.converged);
            NOVA_REQUIRE(ic_solve.message.find("IC(0)") != std::string::npos);
            NOVA_REQUIRE(jacobi.linear_solvers[step].converged);
            NOVA_REQUIRE_NEAR(ic.diagnostics[step].time,
                              jacobi.diagnostics[step].time, 1e-12);
            NOVA_REQUIRE_NEAR(ic.diagnostics[step].heat_content,
                              jacobi.diagnostics[step].heat_content, 2e-9);
            NOVA_REQUIRE_NEAR(ic.diagnostics[step].balance_residual, 0.0, 2e-10);
            NOVA_REQUIRE_NEAR(jacobi.diagnostics[step].balance_residual, 0.0, 2e-10);
        }
    }
}


NOVA_TEST(cached_cg_reuses_identical_matrices_and_rebuilds_changed_coefficients_or_patterns) {
    const auto matrix = CsrMatrix::from_triplets(3, 3, {
        {0,0,4}, {1,1,5}, {2,2,6},
        {0,1,-1}, {1,0,-1}, {1,2,-1}, {2,1,-1}});
    const auto changed_coefficients = CsrMatrix::from_triplets(3, 3, {
        {0,0,4}, {1,1,5.5}, {2,2,6},
        {0,1,-1}, {1,0,-1}, {1,2,-1}, {2,1,-1}});
    // Same dimensions and nnz as matrix, but a different sparsity pattern.
    const auto changed_pattern = CsrMatrix::from_triplets(3, 3, {
        {0,0,4}, {1,1,5}, {2,2,6},
        {0,2,-1}, {2,0,-1}, {1,2,-1}, {2,1,-1}});
    const std::array solution_one{1.0, -2.0, 3.0};
    const std::array solution_two{-4.0, 3.0, 1.0};
    const auto options = ic_options();
    CachedConjugateGradient cache;

    const auto first = cache.solve(matrix, matrix.multiply(solution_one), options);
    NOVA_REQUIRE(first.converged && cache.factor_builds() == 1);
    const auto second = cache.solve(matrix, matrix.multiply(solution_two), options);
    NOVA_REQUIRE(second.converged && cache.factor_builds() == 1);
    for (std::size_t i = 0; i < 3; ++i) {
        NOVA_REQUIRE_NEAR(first.solution[i], solution_one[i], 1e-11);
        NOVA_REQUIRE_NEAR(second.solution[i], solution_two[i], 1e-11);
    }

    const auto third = cache.solve(
        changed_coefficients, changed_coefficients.multiply(solution_one), options);
    NOVA_REQUIRE(third.converged && cache.factor_builds() == 2);
    const auto fourth = cache.solve(
        changed_pattern, changed_pattern.multiply(solution_two), options);
    NOVA_REQUIRE(fourth.converged && cache.factor_builds() == 3);
    for (std::size_t i = 0; i < 3; ++i) {
        NOVA_REQUIRE_NEAR(third.solution[i], solution_one[i], 1e-11);
        NOVA_REQUIRE_NEAR(fourth.solution[i], solution_two[i], 1e-11);
    }

    // A different preconditioner bypasses the IC(0) cache without invalidating it.
    auto jacobi = options;
    jacobi.preconditioner = PreconditionerType::Jacobi;
    NOVA_REQUIRE(cache.solve(changed_pattern, changed_pattern.multiply(solution_one), jacobi).converged);
    NOVA_REQUIRE(cache.factor_builds() == 3);
    NOVA_REQUIRE(cache.solve(changed_pattern, changed_pattern.multiply(solution_one), options).converged);
    NOVA_REQUIRE(cache.factor_builds() == 3);

    // A failed rebuild must leave the last valid matrix/factor pair usable.
    const auto bad_spd = CsrMatrix::from_triplets(4,4,{
        {0,0,4}, {1,1,9}, {2,2,9}, {3,3,9}, {0,1,3}, {1,0,3},
        {1,2,-6}, {2,1,-6}, {2,3,-4}, {3,2,-4}, {0,3,-5}, {3,0,-5}});
    NOVA_REQUIRE(rejects([&] {
        (void)cache.solve(bad_spd, std::array{1.0,1.0,1.0,1.0}, options);
    }));
    NOVA_REQUIRE(cache.factor_builds() == 3);
    NOVA_REQUIRE(cache.solve(changed_pattern, changed_pattern.multiply(solution_two), options).converged);
    NOVA_REQUIRE(cache.factor_builds() == 3);
}
