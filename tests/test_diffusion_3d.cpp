#include "test_framework.hpp"

#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using nova::fem::BoundaryFace3D;
using nova::fem::ScalarDiffusion3DProblem;
using nova::fem::StructuredBoxMesh3D;
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

void append_box_component(
    TetrahedralMesh3D& mesh,
    std::vector<BoundaryFace3D>& boundary_faces,
    const StructuredBoxMesh3D& box,
    const double x_offset
) {
    const std::size_t node_offset = mesh.nodes.size();
    for (const auto& point : box.mesh.nodes) {
        mesh.nodes.push_back({
            point.x + x_offset, point.y, point.z
        });
    }
    for (auto cell : box.mesh.elements) {
        for (std::size_t& node : cell) {
            node += node_offset;
        }
        mesh.elements.push_back(cell);
    }
    for (auto face : box.boundary_faces) {
        for (std::size_t& node : face) {
            node += node_offset;
        }
        boundary_faces.push_back(face);
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

NOVA_TEST(diffusion_3d_solves_balanced_pure_neumann_components_with_zero_mean) {
    auto first = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        2U, 1U, 1U
    );
    for (auto& point : first.mesh.nodes) {
        if (point.x == 0.5) {
            point.x = 0.1;
        }
    }
    const auto second = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        1U, 1U, 1U
    );

    ScalarDiffusion3DProblem problem;
    std::vector<BoundaryFace3D> boundary_faces;
    append_box_component(problem.mesh, boundary_faces, first, 0.0);
    const std::size_t first_component_node_count = problem.mesh.nodes.size();
    const std::size_t first_component_element_count =
        problem.mesh.elements.size();
    append_box_component(problem.mesh, boundary_faces, second, 2.0);
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };

    for (const BoundaryFace3D& face : boundary_faces) {
        const bool on_x_min =
            face_on_x(problem.mesh, face, 0.0) ||
            face_on_x(problem.mesh, face, 2.0);
        const bool on_x_max =
            face_on_x(problem.mesh, face, 1.0) ||
            face_on_x(problem.mesh, face, 3.0);
        const double flux = on_x_min ? 1.0 : (on_x_max ? -1.0 : 0.0);
        problem.neumann_faces.push_back({
            face,
            [flux](double, double, double) { return flux; }
        });
    }

    const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        const double local_x = node < first_component_node_count
            ? problem.mesh.nodes[node].x
            : problem.mesh.nodes[node].x - 2.0;
        NOVA_REQUIRE_NEAR(result.values[node], local_x - 0.5, 2.0e-9);
    }

    const auto quality = nova::fem::evaluate_mesh_quality_3d(problem.mesh);
    long double first_integral = 0.0L;
    long double first_volume = 0.0L;
    long double second_integral = 0.0L;
    long double second_volume = 0.0L;
    for (std::size_t element = 0;
         element < problem.mesh.elements.size();
         ++element) {
        const auto& cell = problem.mesh.elements[element];
        const long double volume = quality.elements[element].volume;
        const long double nodal_sum =
            static_cast<long double>(result.values[cell[0]]) +
            static_cast<long double>(result.values[cell[1]]) +
            static_cast<long double>(result.values[cell[2]]) +
            static_cast<long double>(result.values[cell[3]]);
        if (element < first_component_element_count) {
            first_integral += volume * nodal_sum / 4.0L;
            first_volume += volume;
        } else {
            second_integral += volume * nodal_sum / 4.0L;
            second_volume += volume;
        }
    }
    NOVA_REQUIRE_NEAR(
        static_cast<double>(first_integral / first_volume), 0.0, 2.0e-10
    );
    NOVA_REQUIRE_NEAR(
        static_cast<double>(second_integral / second_volume), 0.0, 2.0e-10
    );
}

NOVA_TEST(diffusion_3d_rejects_incompatible_pure_neumann_loads_per_component) {
    const auto unit_box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        1U, 1U, 1U
    );

    ScalarDiffusion3DProblem unbalanced;
    unbalanced.mesh = unit_box.mesh;
    unbalanced.diffusivity = [](double, double, double) { return 1.0; };
    unbalanced.source = [](double, double, double) { return 1.0; };
    for (const BoundaryFace3D& face : unit_box.boundary_faces) {
        unbalanced.robin_faces.push_back({
            face,
            [](double, double, double) { return 0.0; },
            [](double, double, double) { return 0.0; }
        });
    }

    bool rejected_single_component = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_3d(unbalanced));
    } catch (const std::invalid_argument&) {
        rejected_single_component = true;
    }
    NOVA_REQUIRE(rejected_single_component);

    ScalarDiffusion3DProblem disconnected;
    std::vector<BoundaryFace3D> unused_faces;
    append_box_component(disconnected.mesh, unused_faces, unit_box, 0.0);
    append_box_component(disconnected.mesh, unused_faces, unit_box, 2.0);
    disconnected.diffusivity = [](double, double, double) { return 1.0; };
    disconnected.source = [](const double x, double, double) {
        return x < 1.5 ? 1.0 : -1.0;
    };

    bool rejected_opposite_component_loads = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_3d(disconnected));
    } catch (const std::invalid_argument&) {
        rejected_opposite_component_loads = true;
    }
    NOVA_REQUIRE(rejected_opposite_component_loads);
}

NOVA_TEST(diffusion_3d_positive_robin_transfer_anchors_without_dirichlet_data) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        1U, 1U, 1U
    );

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    for (const BoundaryFace3D& face : box.boundary_faces) {
        problem.robin_faces.push_back({
            face,
            [](double, double, double) { return 2.0; },
            [](double, double, double) { return 3.0; }
        });
    }

    const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (const double value : result.values) {
        NOVA_REQUIRE_NEAR(value, 3.0, 2.0e-10);
    }
}

NOVA_TEST(diffusion_3d_rejects_negative_robin_data) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        1U, 1U, 1U
    );

    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    problem.dirichlet_constraints.push_back({0U, 0.0});
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

NOVA_TEST(diffusion_3d_pure_neumann_residual_describes_the_gauged_cg_iterate) {
    ScalarDiffusion3DProblem problem;
    problem.mesh = {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {{0, 1, 2, 3}}};
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    problem.neumann_faces = {
        {{0, 2, 3}, [](double, double, double) { return 1.0; }},
        {{1, 2, 3}, [](double, double, double) { return -1.0 / std::sqrt(3.0); }}
    };
    const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    const auto& x = result.linear_solver.solution;
    // Pin node zero. The free 3x3 stiffness block is I/6, with load (1,0,0)/6.
    const double residual = std::hypot(std::hypot(x[0], (1.0 - x[1]) / 6.0),
                                       std::hypot(x[2] / 6.0, x[3] / 6.0));
    NOVA_REQUIRE_NEAR(result.linear_solver.residual_norm, residual, 1.0e-14);
    for (std::size_t i = 0; i < x.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.values[i], problem.mesh.nodes[i].x - 0.25, 1.0e-13);
    }
}


// Independent 1D transmission solution embedded in a fitted 3D tetrahedral
// grid. For k=2 on x<1/2 and k=4 on x>1/2, u(0)=0, outward
// q_n(1)=-1, and insulated side walls, the exact solution is:
//   u(x)=x/2                     (x<=1/2)
//   u(x)=1/4+(x-1/2)/4           (x>=1/2)
// Both layers carry the SAME physical flux q_x=-k*du/dx=-1.
// The right-plate potential is 3/8 and the stored energy is 3/16.
// This test validates the continuum interface law independently of
// the FEM assembly and flux postprocessor.
NOVA_TEST(diffusion_3d_two_layer_exact_transmission_and_energy) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    const auto exact = [](const double x) {
        return x <= 0.5 ? 0.5 * x : 0.25 + 0.25 * (x - 0.5);
    };

    // Check both independent material APIs and both representations of
    // insulation (omitted faces versus explicit zero Neumann data).
    for (int representation = 0; representation < 2; ++representation) {
        ScalarDiffusion3DProblem problem;
        problem.mesh = box.mesh;
        if (representation == 0) {
            for (const auto& cell : problem.mesh.elements) {
                double centroid_x = 0.0;
                for (const auto node : cell) {
                    centroid_x += problem.mesh.nodes[node].x / 4.0;
                }
                problem.cell_diffusivity.push_back(
                    centroid_x < 0.5 ? 2.0 : 4.0);
            }
        } else {
            problem.diffusivity = [](const double x, double, double) {
                return x < 0.5 ? 2.0 : 4.0;
            };
        }

        for (const BoundaryFace3D& face : box.boundary_faces) {
            if (face_on_x(problem.mesh, face, 0.0)) {
                problem.dirichlet_faces.push_back({
                    face, [](double, double, double) { return 0.0; }});
            } else if (face_on_x(problem.mesh, face, 1.0)) {
                problem.neumann_faces.push_back({
                    face, [](double, double, double) { return -1.0; }});
            } else if (representation == 0) {
                problem.neumann_faces.push_back({
                    face, [](double, double, double) { return 0.0; }});
            }
        }

        const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
        NOVA_REQUIRE(result.linear_solver.converged);
        for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
            const double x = problem.mesh.nodes[i].x;
            NOVA_REQUIRE_NEAR(result.values[i], exact(x), 1.0e-8);
        }

        const auto fluxes = nova::fem::evaluate_p1_tetrahedral_flux_3d(
            problem, result.values);
        const auto quality = nova::fem::evaluate_mesh_quality_3d(problem.mesh);
        NOVA_REQUIRE(fluxes.size() == problem.mesh.elements.size());
        NOVA_REQUIRE_NEAR(quality.total_volume, 1.0, 1.0e-12);

        double stored_energy = 0.0;
        for (std::size_t i = 0; i < problem.mesh.elements.size(); ++i) {
            const auto& cell = problem.mesh.elements[i];
            double centroid_x = 0.0;
            for (const auto node : cell) {
                centroid_x += problem.mesh.nodes[node].x / 4.0;
            }
            const double k = centroid_x < 0.5 ? 2.0 : 4.0;
            const double analytic_gradient = 1.0 / k;
            NOVA_REQUIRE_NEAR(fluxes[i].gradient.x, analytic_gradient, 1.0e-8);
            NOVA_REQUIRE_NEAR(fluxes[i].gradient.y, 0.0, 1.0e-8);
            NOVA_REQUIRE_NEAR(fluxes[i].gradient.z, 0.0, 1.0e-8);
            NOVA_REQUIRE_NEAR(fluxes[i].flux.x, -1.0, 1.0e-8);
            NOVA_REQUIRE_NEAR(fluxes[i].flux.y, 0.0, 1.0e-8);
            NOVA_REQUIRE_NEAR(fluxes[i].flux.z, 0.0, 1.0e-8);
            const auto& grad = fluxes[i].gradient;
            stored_energy += 0.5 * k *
                (grad.x * grad.x + grad.y * grad.y + grad.z * grad.z) *
                quality.elements[i].volume;
        }
        NOVA_REQUIRE_NEAR(stored_energy, 3.0 / 16.0, 1.0e-8);
        NOVA_REQUIRE_NEAR(2.0 * stored_energy, exact(1.0), 1.0e-8);
    }
}


// Separate source-loaded interface benchmark, not the exact-P1 patch test:
// for k=2 on the left, k=4 on the right, -div(k grad u)=4,
// u(0)=u(1)=0, and insulated sidewalls, flux continuity fixes:
//   u_left(x)  = 5*x/6 - x*x
//   u_right(x) = 1/6 + 5*(x-1/2)/12 - (x*x-1/4)/2.
// The H1 energy-seminorm error is evaluated by exact tetrahedron moments,
// independently of NOVA's residual estimator or numerical quadrature.
NOVA_TEST(diffusion_3d_two_layer_manufactured_source_energy_convergence) {
    constexpr std::array<std::size_t, 3> resolutions{2U, 4U, 8U};
    std::array<double, 3> energy_error_squared{};

    for (std::size_t level = 0; level < resolutions.size(); ++level) {
        const auto n = resolutions[level];
        const auto box = nova::fem::make_box_mesh_3d(
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, n, n);
        ScalarDiffusion3DProblem problem;
        problem.mesh = box.mesh;
        problem.source = [](double, double, double) { return 4.0; };

        for (const auto& cell : problem.mesh.elements) {
            double centroid_x = 0.0;
            for (const auto node : cell) {
                centroid_x += problem.mesh.nodes[node].x / 4.0;
            }
            problem.cell_diffusivity.push_back(
                centroid_x < 0.5 ? 2.0 : 4.0);
        }
        for (const BoundaryFace3D& face : box.boundary_faces) {
            if (face_on_x(problem.mesh, face, 0.0) ||
                face_on_x(problem.mesh, face, 1.0)) {
                problem.dirichlet_faces.push_back({
                    face, [](double, double, double) { return 0.0; }});
            }
        }

        const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
        NOVA_REQUIRE(result.linear_solver.converged);
        const auto fluxes = nova::fem::evaluate_p1_tetrahedral_flux_3d(
            problem, result.values);
        const auto quality = nova::fem::evaluate_mesh_quality_3d(problem.mesh);

        double error_squared = 0.0;
        for (std::size_t i = 0; i < problem.mesh.elements.size(); ++i) {
            double sum_x = 0.0;
            double sum_x_squared = 0.0;
            for (const auto node : problem.mesh.elements[i]) {
                const double x = problem.mesh.nodes[node].x;
                sum_x += x;
                sum_x_squared += x * x;
            }
            const double mean_x = sum_x / 4.0;
            const double mean_x_squared =
                (sum_x * sum_x + sum_x_squared) / 20.0;
            const bool left = mean_x < 0.5;
            const double a = left ? -2.0 : -1.0;
            const double b = left ? 5.0 / 6.0 : 5.0 / 12.0;
            const double k = left ? 2.0 : 4.0;
            const auto& gradient = fluxes[i].gradient;
            const double delta = b - gradient.x;
            const double mean_squared_gradient_error =
                delta * delta + 2.0 * a * delta * mean_x +
                a * a * mean_x_squared +
                gradient.y * gradient.y + gradient.z * gradient.z;
            error_squared += k * quality.elements[i].volume *
                mean_squared_gradient_error;
        }
        NOVA_REQUIRE(std::isfinite(error_squared));
        NOVA_REQUIRE(error_squared > 0.0);
        NOVA_REQUIRE_NEAR(quality.total_volume, 1.0, 1.0e-12);
        energy_error_squared[level] = error_squared;
    }

    // Genuine continuum convergence, not just a shrinking self-estimate.
    NOVA_REQUIRE(energy_error_squared[0] > energy_error_squared[1]);
    NOVA_REQUIRE(energy_error_squared[1] > energy_error_squared[2]);
    NOVA_REQUIRE(energy_error_squared[0] / energy_error_squared[1] > 1.5);
    NOVA_REQUIRE(energy_error_squared[1] / energy_error_squared[2] > 1.5);
}
