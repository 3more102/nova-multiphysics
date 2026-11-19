#include "test_framework.hpp"

#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/error_norms_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using namespace nova::fem;

const ExactScalarField3D quadratic{
    [](double x, double y, double z) { return x*x + y*y + z*z; },
    [](double x, double, double) { return 2.0*x; },
    [](double, double y, double) { return 2.0*y; },
    [](double, double, double z) { return 2.0*z; }
};

std::vector<double> sample_quadratic(const TetrahedralMesh3D& mesh) {
    std::vector<double> u;
    u.reserve(mesh.nodes.size());
    for (const auto& p : mesh.nodes) u.push_back(quadratic.value(p.x, p.y, p.z));
    return u;
}

ErrorNorms3D interpolation_error(std::size_t n) {
    const auto box = make_box_mesh_3d(0,1,0,1,0,1,n,n,n);
    const auto values = sample_quadratic(box.mesh);
    return compute_error_norms_3d(box.mesh, values, quadratic);
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

NOVA_TEST(error_norms_3d_reproduce_affine_p1_exactly_for_both_orientations) {
    TetrahedralMesh3D mesh{
        {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}},
        {{0,1,2,3}}
    };
    const ExactScalarField3D exact{
        [](double x,double y,double z) { return 2+3*x-4*y+5*z; },
        [](double,double,double) { return 3; },
        [](double,double,double) { return -4; },
        [](double,double,double) { return 5; }
    };
    std::vector<double> values;
    for (const auto& p : mesh.nodes) values.push_back(exact.value(p.x,p.y,p.z));
    for (int orientation = 0; orientation < 2; ++orientation) {
        const auto norms = compute_error_norms_3d(mesh, values, exact);
        NOVA_REQUIRE(norms.l2 < 1e-13);
        NOVA_REQUIRE(norms.h1_seminorm < 1e-13);
        NOVA_REQUIRE(norms.h1 < 1e-13);
        std::swap(mesh.elements[0][1], mesh.elements[0][3]);
    }
}

NOVA_TEST(error_norms_3d_quadratic_interpolation_recovers_optimal_spatial_rates) {
    const auto coarse = interpolation_error(2);
    const auto medium = interpolation_error(4);
    const auto fine = interpolation_error(8);
    NOVA_REQUIRE(coarse.l2 > medium.l2 && medium.l2 > fine.l2);
    NOVA_REQUIRE(coarse.h1_seminorm > medium.h1_seminorm);
    NOVA_REQUIRE(medium.h1_seminorm > fine.h1_seminorm);
    // For smooth exact fields, P1 interpolation L2 is O(h^2)
    // and the H1 seminorm is O(h).
    NOVA_REQUIRE(coarse.l2 / medium.l2 > 3.5);
    NOVA_REQUIRE(medium.l2 / fine.l2 > 3.5);
    NOVA_REQUIRE(coarse.h1_seminorm / medium.h1_seminorm > 1.8);
    NOVA_REQUIRE(medium.h1_seminorm / fine.h1_seminorm > 1.8);
}

NOVA_TEST(error_norms_3d_manufactured_poisson_solution_converges) {
    auto solve_measure = [](std::size_t n) {
        const auto box = make_box_mesh_3d(0,1,0,1,0,1,n,n,n);
        ScalarDiffusion3DProblem problem;
        problem.mesh = box.mesh;
        problem.diffusivity = [](double,double,double) { return 1.0; };
        problem.source = [](double,double,double) { return -6.0; };
        for (const auto& face : box.boundary_faces) {
            problem.dirichlet_faces.push_back({face, quadratic.value});
        }
        nova::linalg::SolverOptions options;
        options.relative_tolerance = 1e-12;
        options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
        const auto solved = solve_scalar_diffusion_3d(problem, options);
        NOVA_REQUIRE(solved.linear_solver.converged);
        return compute_error_norms_3d(problem.mesh, solved.values, quadratic);
    };
    const auto coarse = solve_measure(2);
    const auto medium = solve_measure(4);
    const auto fine = solve_measure(8);
    NOVA_REQUIRE(coarse.l2 > medium.l2 && medium.l2 > fine.l2);
    NOVA_REQUIRE(coarse.h1_seminorm > medium.h1_seminorm);
    NOVA_REQUIRE(medium.h1_seminorm > fine.h1_seminorm);
    NOVA_REQUIRE(coarse.l2 / medium.l2 > 2.5);
    NOVA_REQUIRE(medium.l2 / fine.l2 > 2.5);
    NOVA_REQUIRE(coarse.h1_seminorm / medium.h1_seminorm > 1.6);
    NOVA_REQUIRE(medium.h1_seminorm / fine.h1_seminorm > 1.6);
}

NOVA_TEST(error_norms_3d_validate_exact_fields_and_node_data) {
    TetrahedralMesh3D mesh{
        {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}},
        {{0,1,2,3}}
    };
    auto values = sample_quadratic(mesh);
    NOVA_REQUIRE(rejects([&] {
        (void)compute_error_norms_3d(mesh, std::vector<double>{0,1}, quadratic);
    }));
    values[0] = std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects([&] {
        (void)compute_error_norms_3d(mesh, values, quadratic);
    }));
    values[0] = 0.0;
    auto missing = quadratic;
    missing.derivative_z = {};
    NOVA_REQUIRE(rejects([&] {
        (void)compute_error_norms_3d(mesh, values, missing);
    }));
    auto nonfinite = quadratic;
    nonfinite.derivative_x = [](double,double,double) {
        return std::numeric_limits<double>::infinity();
    };
    NOVA_REQUIRE(rejects([&] {
        (void)compute_error_norms_3d(mesh, values, nonfinite);
    }));
}
