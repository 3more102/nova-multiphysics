#include "test_framework.hpp"

#include "nova/fem/p1_flux_3d.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

nova::fem::TetrahedralMesh3D unit_tetrahedron() {
    return {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
}

template<class Action>
bool rejects(Action&& action) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

} // namespace

NOVA_TEST(p1_3d_flux_reproduces_affine_scalar_gradient_in_both_orientations) {
    using namespace nova::fem;
    ScalarDiffusion3DProblem problem;
    problem.mesh = unit_tetrahedron();
    problem.diffusivity = [](double, double, double) { return 2.0; };
    const std::array nodal{1.0, 2.0, 3.0, -2.0}; // u = 1+x+2y-3z
    for (int orientation = 0; orientation < 2; ++orientation) {
        const auto fields = evaluate_p1_tetrahedral_flux_3d(problem, nodal);
        NOVA_REQUIRE(fields.size() == 1U);
        NOVA_REQUIRE_NEAR(fields[0].gradient.x, 1.0, 1e-14);
        NOVA_REQUIRE_NEAR(fields[0].gradient.y, 2.0, 1e-14);
        NOVA_REQUIRE_NEAR(fields[0].gradient.z, -3.0, 1e-14);
        NOVA_REQUIRE_NEAR(fields[0].flux.x, -2.0, 1e-14);
        NOVA_REQUIRE_NEAR(fields[0].flux.y, -4.0, 1e-14);
        NOVA_REQUIRE_NEAR(fields[0].flux.z, 6.0, 1e-14);
        std::swap(problem.mesh.elements[0][0], problem.mesh.elements[0][1]);
    }
}

NOVA_TEST(p1_3d_flux_uses_degree_two_quadrature_for_variable_scalar_material) {
    using namespace nova::fem;
    ScalarDiffusion3DProblem problem;
    problem.mesh = unit_tetrahedron();
    const std::array nodal{1.0, 2.0, 3.0, -2.0};
    problem.diffusivity = [](double x, double, double) { return 1.0 + x; };
    const auto linear = evaluate_p1_tetrahedral_flux_3d(problem, nodal);
    // The mean x-coordinate of this tetrahedron is 1/4.
    NOVA_REQUIRE_NEAR(linear[0].flux.x, -1.25, 1e-14);
    NOVA_REQUIRE_NEAR(linear[0].flux.y, -2.50, 1e-14);
    NOVA_REQUIRE_NEAR(linear[0].flux.z, 3.75, 1e-14);

    // The exact mean of x^2 over the reference tetrahedron is 1/10.
    // A single centroid sample would return 1 + 1/16 instead of 1.1.
    problem.diffusivity = [](double x, double, double) { return 1.0 + x*x; };
    const auto quadratic = evaluate_p1_tetrahedral_flux_3d(problem, nodal);
    NOVA_REQUIRE_NEAR(quadratic[0].flux.x, -1.10, 1e-14);
    NOVA_REQUIRE_NEAR(quadratic[0].flux.y, -2.20, 1e-14);
    NOVA_REQUIRE_NEAR(quadratic[0].flux.z, 3.30, 1e-14);
}

NOVA_TEST(p1_3d_flux_handles_all_six_anisotropic_tensor_components) {
    using namespace nova::fem;
    ScalarDiffusion3DProblem problem;
    problem.mesh = unit_tetrahedron();
    problem.diffusivity_tensor = [](double, double, double) {
        return SymmetricTensor3D{4.0,3.0,2.0,1.0,-0.5,0.75};
    };
    const std::array nodal{1.0,2.0,3.0,0.0}; // u = 1+x+2y-z
    const auto fields = evaluate_p1_tetrahedral_flux_3d(problem, nodal);
    NOVA_REQUIRE_NEAR(fields[0].flux.x, -6.50, 1e-14);
    NOVA_REQUIRE_NEAR(fields[0].flux.y, -6.25, 1e-14);
    NOVA_REQUIRE_NEAR(fields[0].flux.z, 1.00, 1e-14);
}

NOVA_TEST(p1_3d_flux_rejects_invalid_materials_and_nodal_data) {
    using namespace nova::fem;
    ScalarDiffusion3DProblem problem;
    problem.mesh = unit_tetrahedron();
    const std::array nodal{0.0, 1.0, 2.0, 3.0};
    NOVA_REQUIRE(rejects([&] { (void)evaluate_p1_tetrahedral_flux_3d(problem, nodal); }));
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.diffusivity_tensor = [](double, double, double) {
        return SymmetricTensor3D{1,1,1,0,0,0};
    };
    NOVA_REQUIRE(rejects([&] { (void)evaluate_p1_tetrahedral_flux_3d(problem, nodal); }));
    problem.diffusivity = {};
    problem.diffusivity_tensor = [](double, double, double) {
        return SymmetricTensor3D{1,1,1,2,0,0};
    };
    NOVA_REQUIRE(rejects([&] { (void)evaluate_p1_tetrahedral_flux_3d(problem, nodal); }));
    problem.diffusivity_tensor = {};
    problem.diffusivity = [](double, double, double) { return -1.0; };
    NOVA_REQUIRE(rejects([&] { (void)evaluate_p1_tetrahedral_flux_3d(problem, nodal); }));
    problem.diffusivity = [](double, double, double) { return 1.0; };
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_p1_tetrahedral_flux_3d(problem, std::array{1.0,2.0,3.0});
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_p1_tetrahedral_flux_3d(problem,
            std::array{1.0,2.0,3.0,std::numeric_limits<double>::quiet_NaN()});
    }));
}
