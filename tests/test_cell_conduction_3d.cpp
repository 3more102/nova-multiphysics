#include "test_framework.hpp"
#include "nova/fem/p1_flux_3d.hpp"
#include "../src/fem/p1_integration_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {
using namespace nova::fem;

ScalarDiffusion3DProblem single_cell() {
    ScalarDiffusion3DProblem problem;
    problem.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    problem.source = [](double, double, double) { return 0.0; };
    problem.cell_diffusivity = {6.0};
    return problem;
}

template<class F>
bool rejects(F&& action) {
    try { action(); } catch (const std::invalid_argument&) { return true; }
    return false;
}
} // namespace

NOVA_TEST(cell_conductivity_stiffness_matches_all_entries_in_either_winding) {
    auto problem = single_cell();
    // k*V=1. Rows are dot products of gradients (-1,-1,-1), ex, ey, ez.
    const double expected[4][4] = {{3,-1,-1,-1}, {-1,1,0,0}, {-1,0,1,0}, {-1,0,0,1}};
    const double tensor_expected[4][4] = {
        {11.5,-4.5,-4.75,-2.25}, {-4.5,4,1,-0.5}, {-4.75,1,3,0.75}, {-2.25,-0.5,0.75,2}
    };
    for (bool tensor : {false, true}) {
        if (tensor) {
            problem.cell_diffusivity.clear();
            problem.diffusivity_tensor = [](double, double, double) {
                return SymmetricTensor3D{4,3,2,1,-0.5,0.75};
            };
        }
        const detail::DiffusionMaterials3D materials(problem);
        auto& cell = problem.mesh.elements.front();
        for (int orientation = 0; orientation < 2; ++orientation) {
            const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(problem.mesh, 0);
            const auto matrix = detail::integrate_p1_volume_3d(
                problem.mesh, cell, basis, materials.on_cell(0), [](auto&&...) {});
            for (std::size_t i = 0; i < 4; ++i) {
                double row_sum = 0.0;
                for (std::size_t j = 0; j < 4; ++j) {
                    const double exact = tensor ? tensor_expected[cell[i]][cell[j]] / 6.0
                                                : expected[cell[i]][cell[j]];
                    NOVA_REQUIRE_NEAR(matrix[i][j], exact, 2e-15);
                    NOVA_REQUIRE(matrix[i][j] == matrix[j][i]);
                    row_sum += matrix[i][j];
                }
                NOVA_REQUIRE_NEAR(row_sum, 0.0, 2e-15);
            }
            std::swap(cell[0], cell[2]);
        }
    }
}

NOVA_TEST(cell_conductivity_solve_and_export_reject_bad_or_ambiguous_materials) {
    auto problem = single_cell();
    const std::array nodal{0.0, 1.0, 2.0, 3.0};
    auto rejected_by_both = [&] {
        NOVA_REQUIRE(rejects([&] { (void)solve_scalar_diffusion_3d(problem); }));
        NOVA_REQUIRE(rejects([&] { (void)evaluate_p1_tetrahedral_flux_3d(problem, nodal); }));
    };
    for (const double k : {0.0, -2.0, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
        problem.cell_diffusivity = {k};
        rejected_by_both();
    }
    problem.cell_diffusivity = {1.0, 2.0};
    rejected_by_both();
    problem.cell_diffusivity = {1.0};
    problem.diffusivity = [](double, double, double) { return 2.0; };
    rejected_by_both();
    problem.diffusivity = {};
    problem.diffusivity_tensor = [](double, double, double) { return SymmetricTensor3D{1,1,1}; };
    rejected_by_both();
}

NOVA_TEST(cell_flux_preserves_zero_gradient_and_resolvable_extreme_temperature_differences) {
    auto problem = single_cell();
    problem.mesh.nodes[2] = {0.27, 1.13, 0.39};
    problem.mesh.nodes[3] = {0.17, 0.37, 0.91};
    const std::array constant{1e280, 1e280, 1e280, 1e280};
    const auto uniform = evaluate_p1_tetrahedral_flux_3d(problem, constant).front().flux;
    NOVA_REQUIRE(uniform.x == 0 && uniform.y == 0 && uniform.z == 0);

    problem = single_cell();
    problem.mesh.nodes[1].x = 4.0;
    problem.cell_diffusivity = {0.5};
    const std::array hot_cold{-1e308, 1e308, -1e308, -1e308};
    const auto field = evaluate_p1_tetrahedral_flux_3d(problem, hot_cold).front();
    NOVA_REQUIRE_NEAR(field.gradient.x / 5e307, 1.0, 1e-15);
    NOVA_REQUIRE_NEAR(field.flux.x / -2.5e307, 1.0, 1e-15);
}

NOVA_TEST(bonded_oblique_interface_transmits_normal_flux_and_preserves_tangential_jump) {
    auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 4, 4, 4);
    ScalarDiffusion3DProblem problem;
    problem.mesh = std::move(box.mesh);
    for (auto& cell : problem.mesh.elements) {
        double mean_x = 0;
        for (auto node : cell) mean_x += problem.mesh.nodes[node].x / 4;
        problem.cell_diffusivity.push_back(mean_x < 0.5 ? 2.0 : 6.0);
        if (cell[0] % 2 == 0) std::swap(cell[1], cell[3]);
    }
    // Orthogonal rotation with first column (2,-1,2)/3 normal to the interface.
    auto rotate = [](double u, double v, double w) -> Point3D {
        return {(2*u+2*v-w)/3, (-u+2*v+2*w)/3, (2*u-v+2*w)/3};
    };
    for (auto& p : problem.mesh.nodes) p = rotate(p.x, p.y, p.z);
    auto exact = [](double x, double y, double z) {
        const double u = (2*x-y+2*z)/3, v = (2*x+2*y-z)/3, w = (-x+2*y+2*z)/3;
        return 300 + 3*v - 4*w - 12*(std::min(u, 0.5)/2 + std::max(u-0.5, 0.0)/6);
    };
    for (const auto& face : box.boundary_faces) problem.dirichlet_faces.push_back({face, exact});
    problem.source = [](double, double, double) { return 0.0; };
    nova::linalg::SolverOptions options;
    options.relative_tolerance = 1e-13;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    const auto result = solve_scalar_diffusion_3d(problem, options);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t i = 0; i < result.values.size(); ++i) {
        const auto& p = problem.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(result.values[i], exact(p.x, p.y, p.z), 2e-9);
    }
    const auto fields = evaluate_p1_tetrahedral_flux_3d(problem, result.values);
    for (std::size_t e = 0; e < fields.size(); ++e) {
        const double k = problem.cell_diffusivity[e];
        const auto expected = rotate(12, -3*k, 4*k);
        NOVA_REQUIRE_NEAR(fields[e].flux.x, expected.x, 2e-8);
        NOVA_REQUIRE_NEAR(fields[e].flux.y, expected.y, 2e-8);
        NOVA_REQUIRE_NEAR(fields[e].flux.z, expected.z, 2e-8);
    }
}


NOVA_TEST(cell_heat_source_adds_to_spatial_density_and_rejects_invalid_data) {
    auto problem = single_cell();
    problem.dirichlet_constraints = {{0, 0.0}};
    problem.source = [](double, double, double) { return 40.0; };
    problem.cell_source = {80.0};
    const auto combined = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(combined.linear_solver.converged);

    problem.source = [](double, double, double) { return 120.0; };
    problem.cell_source.clear();
    const auto spatial = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(spatial.linear_solver.converged);
    for (std::size_t node = 0; node < combined.values.size(); ++node) {
        NOVA_REQUIRE_NEAR(combined.values[node], spatial.values[node], 1e-11);
    }

    problem.source = [](double, double, double) { return 240.0; };
    problem.cell_source = {-120.0};
    const auto cooling = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(cooling.linear_solver.converged);
    for (std::size_t node = 0; node < spatial.values.size(); ++node) {
        NOVA_REQUIRE_NEAR(cooling.values[node], spatial.values[node], 1e-11);
    }

    problem.cell_source = {1.0, 2.0};
    NOVA_REQUIRE(rejects([&] { (void)solve_scalar_diffusion_3d(problem); }));
    for (const double bad : {std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
        problem.cell_source = {bad};
        NOVA_REQUIRE(rejects([&] { (void)solve_scalar_diffusion_3d(problem); }));
    }
}

NOVA_TEST(cell_heat_source_matches_independent_spatial_piecewise_source) {
    const auto box = make_box_mesh_3d(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2, 2, 2);
    ScalarDiffusion3DProblem tagged;
    tagged.mesh = box.mesh;
    tagged.source = [](double, double, double) { return 0.0; };
    for (const auto& cell : tagged.mesh.elements) {
        double center_x = 0.0;
        for (auto node : cell) center_x += tagged.mesh.nodes[node].x / 4.0;
        tagged.cell_diffusivity.push_back(center_x < 0.5 ? 2.0 : 6.0);
        tagged.cell_source.push_back(center_x < 0.5 ? 120.0 : 360.0);
    }
    const ScalarField3D hot = [](double, double, double) { return 400.0; };
    const ScalarField3D h = [](double, double, double) { return 12.0; };
    const ScalarField3D ambient = [](double, double, double) { return 300.0; };
    for (const auto& face : box.boundary_faces) {
        bool at_hot = true, at_cold = true;
        for (auto node : face) {
            at_hot = at_hot && tagged.mesh.nodes[node].x == 0.0;
            at_cold = at_cold && tagged.mesh.nodes[node].x == 1.0;
        }
        if (at_hot) tagged.dirichlet_faces.push_back({face, hot});
        if (at_cold) tagged.robin_faces.push_back({face, h, ambient});
    }
    nova::linalg::SolverOptions options;
    options.relative_tolerance = 1e-12;
    options.absolute_tolerance = 0.0;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    const auto from_tags = solve_scalar_diffusion_3d(tagged, options);
    NOVA_REQUIRE(from_tags.linear_solver.converged);

    auto from_callback = tagged;
    from_callback.cell_source.clear();
    from_callback.source = [](double x, double, double) {
        return x < 0.5 ? 120.0 : 360.0;
    };
    const auto reference = solve_scalar_diffusion_3d(from_callback, options);
    NOVA_REQUIRE(reference.linear_solver.converged);
    for (std::size_t i = 0; i < from_tags.values.size(); ++i) {
        NOVA_REQUIRE_NEAR(from_tags.values[i], reference.values[i], 2e-9);
    }
}

NOVA_TEST(cell_tensor_conductivity_stiffness_is_symmetric_and_orientation_independent) {
    auto problem = single_cell();
    problem.cell_diffusivity.clear();
    problem.cell_diffusivity_tensors = {{4.0, 3.0, 2.0, 1.0, -0.5, 0.75}};
    const double k[4][4] = {
        {11.5,-4.5,-4.75,-2.25},
        {-4.5,4.0,1.0,-0.5},
        {-4.75,1.0,3.0,0.75},
        {-2.25,-0.5,0.75,2.0}
    };
    const detail::DiffusionMaterials3D materials(problem);
    auto& cell = problem.mesh.elements.front();
    for (int orientation = 0; orientation < 2; ++orientation) {
        const auto gradients = evaluate_tetrahedron_p1_physical_gradients_3d(problem.mesh, 0);
        const auto stiffness = detail::integrate_p1_volume_3d(
            problem.mesh, cell, gradients, materials.on_cell(0), [](auto&&...) {});
        for (std::size_t i = 0; i < 4; ++i) {
            double sum = 0;
            for (std::size_t j = 0; j < 4; ++j) {
                NOVA_REQUIRE_NEAR(stiffness[i][j], k[cell[i]][cell[j]] / 6.0, 2e-15);
                NOVA_REQUIRE(stiffness[i][j] == stiffness[j][i]);
                sum += stiffness[i][j];
            }
            NOVA_REQUIRE_NEAR(sum, 0.0, 2e-15);
        }
        std::swap(cell[0], cell[2]);
    }
}

NOVA_TEST(bonded_mixed_anisotropic_layers_preserve_normal_flux_and_transverse_jump) {
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 4, 2, 2);
    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    const SymmetricTensor3D left{2, 3, 4, 0.4, -0.3, 0.2};
    const SymmetricTensor3D right{6, 7, 9, -0.7, 0.6, -0.4};
    for (auto& cell : problem.mesh.elements) {
        double center_x = 0;
        for (auto id : cell) center_x += problem.mesh.nodes[id].x / 4.0;
        problem.cell_diffusivity_tensors.push_back(center_x < 0.5 ? left : right);
        if (cell[0] % 2 == 0) std::swap(cell[1], cell[3]);
    }
    const auto exact = [](double x, double, double) {
        return 400.0 - 240.0*(std::min(x, 0.5)/2.0 + std::max(x-0.5, 0.0)/6.0);
    };
    for (const auto& face : box.boundary_faces) {
        problem.dirichlet_faces.push_back({face, exact});
    }
    problem.source = [](double, double, double) { return 0.0; };
    nova::linalg::SolverOptions options;
    options.relative_tolerance = 1e-13;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    const auto solved = solve_scalar_diffusion_3d(problem, options);
    NOVA_REQUIRE(solved.linear_solver.converged);
    for (std::size_t i = 0; i < solved.values.size(); ++i) {
        const auto& p = problem.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(solved.values[i], exact(p.x, p.y, p.z), 2e-9);
    }
    const auto result = evaluate_p1_tetrahedral_flux_3d(problem, solved.values);
    for (std::size_t i = 0; i < result.size(); ++i) {
        const auto& tensor = problem.cell_diffusivity_tensors[i];
        NOVA_REQUIRE_NEAR(result[i].flux.x, 240.0, 2e-8);
        NOVA_REQUIRE_NEAR(result[i].flux.y, 240.0*tensor.xy/tensor.xx, 2e-8);
        NOVA_REQUIRE_NEAR(result[i].flux.z, 240.0*tensor.xz/tensor.xx, 2e-8);
    }
}

NOVA_TEST(cell_tensor_rejects_non_spd_wrong_size_and_conflicting_material_data) {
    auto problem = single_cell();
    problem.cell_diffusivity.clear();
    const std::array nodal{0.0, 1.0, 2.0, 3.0};
    auto rejected_by_solver_and_postprocessor = [&] {
        NOVA_REQUIRE(rejects([&] { (void)solve_scalar_diffusion_3d(problem); }));
        NOVA_REQUIRE(rejects([&] { (void)evaluate_p1_tetrahedral_flux_3d(problem, nodal); }));
    };
    for (const auto bad : {
            SymmetricTensor3D{1,1,1,2,0,0},
            SymmetricTensor3D{1,1,1,1,1,1},
            SymmetricTensor3D{0,1,1,0,0,0},
            SymmetricTensor3D{1,1,1,std::numeric_limits<double>::quiet_NaN(),0,0},
            SymmetricTensor3D{1,1,std::numeric_limits<double>::infinity(),0,0,0}}) {
        problem.cell_diffusivity_tensors = {bad};
        rejected_by_solver_and_postprocessor();
    }
    problem.cell_diffusivity_tensors = {{1,1,1}, {2,2,2}};
    rejected_by_solver_and_postprocessor();
    problem.cell_diffusivity_tensors = {{1,1,1}};
    problem.cell_diffusivity = {1};
    rejected_by_solver_and_postprocessor();
    problem.cell_diffusivity.clear();
    problem.diffusivity = [](double,double,double) { return 1.0; };
    rejected_by_solver_and_postprocessor();
    problem.diffusivity = {};
    problem.diffusivity_tensor = [](double,double,double) { return SymmetricTensor3D{1,1,1}; };
    rejected_by_solver_and_postprocessor();
}
