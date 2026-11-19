#include "test_framework.hpp"

#include "nova/fem/adaptive_heat_3d.hpp"
#include "nova/fem/diffusion_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace {
using namespace nova::fem;

double one(double, double, double) { return 1.0; }

TetrahedralMesh3D reference_tet() {
    return {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
}

SymmetricTensor3D material(double scale) {
    return {4*scale, 3*scale, 2*scale, scale, -0.5*scale, 0.75*scale};
}

bool on_plane(const TetrahedralMesh3D& mesh, const BoundaryFace3D& face,
              double Point3D::*coordinate, double level) {
    return std::all_of(face.begin(), face.end(), [&](std::size_t i) {
        return mesh.nodes[i].*coordinate == level;
    });
}

template<class F>
bool rejects(F&& action) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}
} // namespace

NOVA_TEST(tensor_stiffness_matches_all_analytical_entries_in_both_orientations) {
    auto mesh = reference_tet();
    // Six K components, global gradients (-1,-1,-1), ex, ey, ez.
    constexpr double expected[4][4] = {
        {11.5,-4.5,-4.75,-2.25}, {-4.5,4,1,-0.5},
        {-4.75,1,3,0.75}, {-2.25,-0.5,0.75,2}
    };
    for (int orientation = 0; orientation < 2; ++orientation) {
        auto operators = assemble_anisotropic_h1_operators_3d(mesh, one,
            [](double x, double y, double z) { return material(1 + x + 2*y + 3*z); });
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = 0; j < 4; ++j) {
                // The mean coefficient multiplier is 2.5, the volume is 1/6.
                NOVA_REQUIRE_NEAR(operators.stiffness.at(i,j), 2.5/6 * expected[i][j], 4e-15);
                NOVA_REQUIRE(operators.stiffness.at(i,j) == operators.stiffness.at(j,i));
                NOVA_REQUIRE_NEAR(operators.mass.at(i,j), (i == j ? 2.0 : 1.0)/120, 1e-16);
            }
        }
        std::swap(mesh.elements[0][0], mesh.elements[0][2]);
    }
}

NOVA_TEST(tensor_assembly_is_invariant_to_material_units_and_resolves_diagonal_contrast) {
    const auto mesh = reference_tet();
    for (double scale : {1e-200, 1.0, 1e200}) {
        const auto tensor = assemble_anisotropic_h1_operators_3d(mesh, one,
            [scale](double, double, double) { return material(scale); });
        NOVA_REQUIRE_NEAR(tensor.stiffness.at(2,3)/scale, 0.75/6, 2e-16);
        NOVA_REQUIRE_NEAR(tensor.stiffness.at(0,0)/scale, 11.5/6, 3e-15);
    }
    const auto contrast = assemble_anisotropic_h1_operators_3d(mesh, one,
        [](double, double, double) { return SymmetricTensor3D{1e-240,1e240,1,0.25,0,0}; });
    NOVA_REQUIRE_NEAR(contrast.stiffness.at(1,1)/1e-240, 1.0/6, 2e-16);
    NOVA_REQUIRE_NEAR(contrast.stiffness.at(2,2)/1e240, 1.0/6, 2e-16);
    NOVA_REQUIRE_NEAR(contrast.stiffness.at(1,2), 0.25/6, 1e-16);
}

NOVA_TEST(isotropic_tensor_matches_the_scalar_thermal_operators) {
    const auto mesh = make_box_mesh_3d(0,1,0,2,-1,1, 2,2,2).mesh;
    const auto k = [](double x, double y, double z) { return 3 + x*x + y + z; };
    const auto scalar = assemble_scalar_h1_operators_3d(mesh, one, k);
    const auto tensor = assemble_anisotropic_h1_operators_3d(mesh, one,
        [&](double x, double y, double z) {
            const double value = k(x,y,z);
            return SymmetricTensor3D{value,value,value,0,0,0};
        });
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
        for (std::size_t j = 0; j < mesh.nodes.size(); ++j) {
            NOVA_REQUIRE_NEAR(tensor.stiffness.at(i,j), scalar.stiffness.at(i,j), 2e-14);
            NOVA_REQUIRE(tensor.mass.at(i,j) == scalar.mass.at(i,j));
        }
    }
}

NOVA_TEST(scalar_3d_stiffness_keeps_representable_entries_on_rescaled_tetrahedra) {
    auto mesh = reference_tet();
    for (auto& p : mesh.nodes) {
        p.x *= 1e-90;
        p.y *= 1e-90;
        p.z *= 1e-90;
    }
    const auto operators = assemble_scalar_h1_operators_3d(mesh, one,
        [](double,double,double) { return 1e-90; });
    // K11 = k h / 6 even though k times the volume underflows.
    NOVA_REQUIRE_NEAR(operators.stiffness.at(1,1)/1e-180, 1.0/6, 2e-16);

    ScalarDiffusion3DProblem steady;
    steady.mesh = mesh;
    steady.diffusivity = [](double,double,double) { return 1e200; };
    steady.source = [](double,double,double) { return 0.0; };
    steady.dirichlet_constraints = {{1,1}, {2,1}, {3,1}};
    nova::linalg::SolverOptions options;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    // k |grad N|^2 overflows before multiplication by h^3 on the parent.
    const auto solution = solve_scalar_diffusion_3d(steady, options);
    NOVA_REQUIRE(solution.linear_solver.converged);
    NOVA_REQUIRE_NEAR(solution.values[0], 1.0, 2e-15);
}

NOVA_TEST(tensor_diffusion_and_heat_recover_mixed_fluxes_on_a_distorted_cube) {
    auto box = make_box_mesh_3d(0,1,0,1,0,1, 2,2,2);
    for (auto& p : box.mesh.nodes) {
        if (p.x == 0.5 && p.y == 0.5 && p.z == 0.5) p = {0.42,0.63,0.55};
    }
    for (std::size_t e = 0; e < box.mesh.elements.size(); e += 3) {
        std::swap(box.mesh.elements[e][0], box.mesh.elements[e][1]);
    }
    const auto exact = [](double x, double y, double z) { return 1 + x + 2*y - z; };
    ScalarDiffusion3DProblem steady;
    steady.mesh = box.mesh;
    steady.diffusivity_tensor = [](double x, double, double) { return material(1+x); };
    steady.source = [](double, double, double) { return -6.5; };
    TransientHeat3DProblem heat;
    heat.mesh = box.mesh;
    heat.conductivity_tensor = steady.diffusivity_tensor;
    heat.volumetric_capacity = [](double, double, double) { return 2.0; };
    heat.source = [](double, double, double, double) { return -4.5; };
    heat.initial_temperature = exact;
    // K grad(u) = (1+x) (6.5, 6.25, -1). Use the analytic outward flux.
    for (const auto& face : box.boundary_faces) {
        if (on_plane(box.mesh, face, &Point3D::x, 0)) {
            steady.dirichlet_faces.push_back({face, exact});
            heat.dirichlet_faces.push_back({face, [exact](double x, double y, double z, double t) {
                return exact(x,y,z) + t;
            }});
        } else if (on_plane(box.mesh, face, &Point3D::x, 1)) {
            const auto ambient = [exact](double x, double y, double z) { return exact(x,y,z) + 13.0/4; };
            steady.robin_faces.push_back({face, [](double,double,double) { return 4.0; }, ambient});
            heat.robin_faces.push_back({face, [](double,double,double,double) { return 4.0; },
                [ambient](double x,double y,double z,double t) { return ambient(x,y,z) + t; }});
        } else {
            double multiplier = 1.0; // z=1
            if (on_plane(box.mesh, face, &Point3D::z, 0)) multiplier = -1.0;
            if (on_plane(box.mesh, face, &Point3D::y, 0)) multiplier = 6.25;
            if (on_plane(box.mesh, face, &Point3D::y, 1)) multiplier = -6.25;
            steady.neumann_faces.push_back({face,
                [multiplier](double x,double,double) { return multiplier*(1+x); }});
            heat.neumann_faces.push_back({face,
                [multiplier](double x,double,double,double) { return multiplier*(1+x); }});
        }
    }
    const auto solution = solve_scalar_diffusion_3d(steady);
    NOVA_REQUIRE(solution.linear_solver.converged);
    for (std::size_t i = 0; i < steady.mesh.nodes.size(); ++i) {
        const auto& p = steady.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(solution.values[i], exact(p.x,p.y,p.z), 2e-9);
    }
    for (auto scheme : {TransientTimeScheme::BackwardEuler, TransientTimeScheme::CrankNicolson}) {
        TransientHeat3DOptions options;
        options.time_step = 0.05;
        options.steps = 3;
        options.scheme = scheme;
        options.solver_options.relative_tolerance = 1e-13;
        const auto result = solve_transient_heat_3d(heat, options);
        for (std::size_t i = 0; i < heat.mesh.nodes.size(); ++i) {
            const auto& p = heat.mesh.nodes[i];
            NOVA_REQUIRE_NEAR(result.temperature[i], exact(p.x,p.y,p.z) + result.final_time, 2e-11);
        }
        for (const auto& d : result.diagnostics) {
            NOVA_REQUIRE_NEAR(d.source_energy, -4.5*options.time_step, 1e-12);
            NOVA_REQUIRE_NEAR(d.neumann_energy, 0, 1e-12);
            NOVA_REQUIRE_NEAR(d.robin_energy, 13*options.time_step, 1e-11);
            NOVA_REQUIRE_NEAR(d.dirichlet_energy, -6.5*options.time_step, 1e-11);
            NOVA_REQUIRE_NEAR(d.heat_content, 4 + 2*d.time, 1e-11);
            NOVA_REQUIRE_NEAR(d.balance_residual, 0, 1e-11);
        }
    }
}

NOVA_TEST(layered_tensor_diffusion_preserves_interface_flux_and_the_floating_mean) {
    const auto box = make_box_mesh_3d(0,1,0,1,0,1, 4,2,2);
    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity_tensor = [](double x, double, double) {
        const double k = x < 0.5 ? 2.0 : 8.0;
        return SymmetricTensor3D{k,4,3,k/4,-k/8,0.5};
    };
    problem.source = [](double, double, double) { return 0.0; };
    for (const auto& face : box.boundary_faces) {
        double flux = 0.125; // z=1
        if (on_plane(box.mesh, face, &Point3D::x, 0)) flux = 1;
        if (on_plane(box.mesh, face, &Point3D::x, 1)) flux = -1;
        if (on_plane(box.mesh, face, &Point3D::y, 0)) flux = 0.25;
        if (on_plane(box.mesh, face, &Point3D::y, 1)) flux = -0.25;
        if (on_plane(box.mesh, face, &Point3D::z, 0)) flux = -0.125;
        problem.neumann_faces.push_back({face, [flux](double,double,double) { return flux; }});
    }
    const auto result = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t i = 0; i < box.mesh.nodes.size(); ++i) {
        const double x = box.mesh.nodes[i].x;
        const double exact = (x <= 0.5 ? x/2 : 0.25 + (x-0.5)/8) - 13.0/64;
        NOVA_REQUIRE_NEAR(result.values[i], exact, 2e-10);
    }
}

NOVA_TEST(conductivity_rejects_invalid_or_ambiguous_materials_before_elimination) {
    ScalarDiffusion3DProblem p;
    p.mesh = reference_tet();
    p.source = one;
    for (std::size_t i = 0; i < 4; ++i) p.dirichlet_constraints.push_back({i, 0});
    NOVA_REQUIRE(rejects([&] { (void)solve_scalar_diffusion_3d(p); }));
    p.diffusivity = one;
    p.diffusivity_tensor = [](double,double,double) { return material(1); };
    NOVA_REQUIRE(rejects([&] { (void)solve_scalar_diffusion_3d(p); }));
    p.diffusivity = {};
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto& invalid : std::array{
             SymmetricTensor3D{0,1,1,0,0,0}, SymmetricTensor3D{1,-1,1,0,0,0},
             SymmetricTensor3D{1,1,nan,0,0,0}, SymmetricTensor3D{1,1,1,0,nan,0},
             SymmetricTensor3D{1,1,1,0,0,2}, SymmetricTensor3D{1,1,1,0.5,0.5,-0.5},
             SymmetricTensor3D{1,1,1,0.7,0.7,-0.7}}) {
        // Only one of the four samples fails; averaging first could hide it.
        p.diffusivity_tensor = [invalid](double x,double,double) { return x > 0.5 ? invalid : material(1); };
        NOVA_REQUIRE(rejects([&] { (void)solve_scalar_diffusion_3d(p); }));
        NOVA_REQUIRE(rejects([&] { (void)assemble_anisotropic_h1_operators_3d(p.mesh, one, p.diffusivity_tensor); }));
    }
    TransientHeat3DProblem heat;
    heat.mesh = p.mesh;
    heat.conductivity = heat.volumetric_capacity = heat.initial_temperature = one;
    heat.conductivity_tensor = p.diffusivity_tensor;
    heat.source = [](double,double,double,double) { return 0.0; };
    NOVA_REQUIRE(rejects([&] { (void)solve_transient_heat_3d(heat); }));
    NOVA_REQUIRE(rejects([&] { (void)solve_adaptive_heat_3d(heat); }));
}
