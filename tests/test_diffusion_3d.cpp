#include "test_framework.hpp"

#include "nova/fem/diffusion_3d.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using nova::fem::BoundaryFace3D;
using nova::fem::ScalarDiffusion3DProblem;
using nova::fem::TetrahedralMesh3D;

[[nodiscard]] bool face_on_x(
    const TetrahedralMesh3D& mesh,
    const BoundaryFace3D& face,
    const double x
) {
    return mesh.nodes[face[0]].x == x &&
        mesh.nodes[face[1]].x == x &&
        mesh.nodes[face[2]].x == x;
}

void prescribe_all_boundary_faces(
    ScalarDiffusion3DProblem& problem,
    const std::function<double(double, double, double)>& exact
) {
    for (const BoundaryFace3D& face :
         nova::fem::extract_boundary_faces_3d(problem.mesh)) {
        problem.dirichlet_faces.push_back({face, exact});
    }
}

[[nodiscard]] double solve_cubic_x_problem(
    const std::size_t subdivisions
) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        subdivisions, subdivisions, subdivisions
    );

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](const double x, double, double) {
        return 1.0 + 3.0 * x;
    };

    for (const BoundaryFace3D& face : box.boundary_faces) {
        if (face_on_x(problem.mesh, face, 0.0) ||
            face_on_x(problem.mesh, face, 1.0)) {
            problem.dirichlet_faces.push_back({
                face,
                [](double, double, double) { return 0.0; }
            });
        }
    }

    const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);

    double maximum_error = 0.0;
    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        const double x = problem.mesh.nodes[node].x;
        const double exact = x - 0.5 * x * x - 0.5 * x * x * x;
        maximum_error = std::max(
            maximum_error, std::abs(result.values[node] - exact)
        );
    }
    return maximum_error;
}

} // namespace

NOVA_TEST(diffusion_3d_reproduces_affine_harmonic_field) {
    const auto box = nova::fem::make_box_mesh_3d(
        -1.0, 2.0, 0.0, 2.0, -0.5, 1.5,
        2U, 2U, 2U
    );

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 2.0; };
    problem.source = [](double, double, double) { return 0.0; };
    const auto exact = [](const double x, const double y, const double z) {
        return 1.5 + 2.0 * x - 3.0 * y + 4.0 * z;
    };
    prescribe_all_boundary_faces(problem, exact);

    const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE(result.values.size() == problem.mesh.nodes.size());

    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        const auto& point = problem.mesh.nodes[node];
        NOVA_REQUIRE_NEAR(
            result.values[node], exact(point.x, point.y, point.z), 2.0e-10
        );
    }
}

NOVA_TEST(diffusion_3d_supports_variable_diffusivity_and_source) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        2U, 2U, 2U
    );

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](const double x, double, double) {
        return 1.0 + x;
    };
    problem.source = [](double, double, double) { return -1.0; };
    prescribe_all_boundary_faces(
        problem,
        [](const double x, double, double) { return x; }
    );

    const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        NOVA_REQUIRE_NEAR(
            result.values[node], problem.mesh.nodes[node].x, 2.0e-10
        );
    }
}

NOVA_TEST(diffusion_3d_solves_mixed_dirichlet_neumann_problem) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        2U, 2U, 2U
    );

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };

    for (const BoundaryFace3D& face : box.boundary_faces) {
        if (face_on_x(problem.mesh, face, 0.0)) {
            problem.dirichlet_faces.push_back({
                face,
                [](double, double, double) { return 0.0; }
            });
        } else {
            const double flux = face_on_x(problem.mesh, face, 1.0) ? -1.0 : 0.0;
            problem.neumann_faces.push_back({
                face,
                [flux](double, double, double) { return flux; }
            });
        }
    }

    const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        NOVA_REQUIRE_NEAR(
            result.values[node], problem.mesh.nodes[node].x, 2.0e-10
        );
    }
}

NOVA_TEST(diffusion_3d_solves_mixed_dirichlet_robin_problem) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        2U, 2U, 2U
    );

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };

    for (const BoundaryFace3D& face : box.boundary_faces) {
        if (face_on_x(problem.mesh, face, 0.0)) {
            problem.dirichlet_faces.push_back({
                face,
                [](double, double, double) { return 0.0; }
            });
        } else if (face_on_x(problem.mesh, face, 1.0)) {
            problem.robin_faces.push_back({
                face,
                [](double, double, double) { return 2.0; },
                [](double, double, double) { return 1.5; }
            });
        } else {
            problem.neumann_faces.push_back({
                face,
                [](double, double, double) { return 0.0; }
            });
        }
    }

    const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        NOVA_REQUIRE_NEAR(
            result.values[node], problem.mesh.nodes[node].x, 2.0e-10
        );
    }
}

NOVA_TEST(diffusion_3d_accepts_reversed_tetrahedron_orientation) {
    auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        2U, 2U, 2U
    );
    std::swap(box.mesh.elements[0][1], box.mesh.elements[0][3]);

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    prescribe_all_boundary_faces(
        problem,
        [](const double x, const double y, const double z) {
            return 1.0 + x - 2.0 * y + 3.0 * z;
        }
    );

    const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        const auto& point = problem.mesh.nodes[node];
        NOVA_REQUIRE_NEAR(
            result.values[node],
            1.0 + point.x - 2.0 * point.y + 3.0 * point.z,
            2.0e-10
        );
    }
}

NOVA_TEST(diffusion_3d_has_second_order_nodal_convergence) {
    const double coarse_error = solve_cubic_x_problem(2U);
    const double medium_error = solve_cubic_x_problem(4U);
    const double fine_error = solve_cubic_x_problem(8U);

    NOVA_REQUIRE(coarse_error > medium_error);
    NOVA_REQUIRE(medium_error > fine_error);
    NOVA_REQUIRE(coarse_error / medium_error > 2.5);
    NOVA_REQUIRE(medium_error / fine_error > 2.5);
}

NOVA_TEST(diffusion_3d_rejects_invalid_inputs_and_boundary_faces) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        1U, 1U, 1U
    );

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 0.0; };
    problem.source = [](double, double, double) { return 0.0; };
    problem.dirichlet_constraints.push_back({0U, 0.0});

    bool rejected_diffusivity = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_3d(problem));
    } catch (const std::invalid_argument&) {
        rejected_diffusivity = true;
    }
    NOVA_REQUIRE(rejected_diffusivity);

    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.neumann_faces.push_back({
        BoundaryFace3D{
            box.mesh.elements[0][0],
            box.mesh.elements[0][1],
            box.mesh.elements[0][3]
        },
        [](double, double, double) { return 0.0; }
    });

    bool rejected_interior_face = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_3d(problem));
    } catch (const std::invalid_argument&) {
        rejected_interior_face = true;
    }
    NOVA_REQUIRE(rejected_interior_face);
}

NOVA_TEST(diffusion_3d_rejects_pure_neumann_and_negative_robin_data) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        1U, 1U, 1U
    );

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    problem.neumann_faces.push_back({
        box.boundary_faces[0],
        [](double, double, double) { return 0.0; }
    });

    bool rejected_pure_neumann = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_3d(problem));
    } catch (const std::invalid_argument&) {
        rejected_pure_neumann = true;
    }
    NOVA_REQUIRE(rejected_pure_neumann);

    problem.dirichlet_constraints.push_back({0U, 0.0});
    problem.neumann_faces.clear();
    problem.robin_faces.push_back({
        box.boundary_faces[0],
        [](double, double, double) { return -1.0; },
        [](double, double, double) { return 0.0; }
    });

    bool rejected_negative_robin = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_3d(problem));
    } catch (const std::invalid_argument&) {
        rejected_negative_robin = true;
    }
    NOVA_REQUIRE(rejected_negative_robin);
}
