#include "test_framework.hpp"

#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/interior_flux_jumps_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace {

// Independent analytic transmission benchmark:
//
// -div(k grad u) = 0 on the unit cube;
// k = 2 on x < 1/2, k = 6 on x > 1/2;
// u = 0 at x = 0; outward q = -1 at x = 1;
// the remaining faces are insulated.
//
// Continuity of u and k*u_x gives u = x/2 on the left and
// u = 1/4 + (x-1/2)/6 on the right. The exact elementwise flux is
// (-1,0,0), the right-end potential is 1/3, and
// the stored energy per unit cross-section is 1/6.
[[nodiscard]] double exact_temperature(const double x) {
    return x <= 0.5 ? 0.5 * x : 0.25 + (x - 0.5) / 6.0;
}

[[nodiscard]] bool all_nodes_on_x(
    const nova::fem::TetrahedralMesh3D& mesh,
    const nova::fem::BoundaryFace3D& face,
    const double x
) {
    for (const auto node : face) {
        if (mesh.nodes[node].x != x) return false;
    }
    return true;
}

} // namespace

NOVA_TEST(diffusion_3d_exact_two_layer_transmission_and_flux_continuity) {
    using namespace nova::fem;

    for (const std::size_t n : std::array<std::size_t, 2>{2U, 4U}) {
        const auto box = make_box_mesh_3d(
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, n, n);
        const auto quality = evaluate_mesh_quality_3d(box.mesh);
        NOVA_REQUIRE_NEAR(quality.total_volume, 1.0, 1.0e-12);

        std::vector<double> exact_cell_k;
        exact_cell_k.reserve(box.mesh.elements.size());
        for (const auto& tetra : box.mesh.elements) {
            double center_x = 0.0;
            for (const auto node : tetra) {
                center_x += box.mesh.nodes[node].x / 4.0;
            }
            // Both grids explicitly fit the dielectric/conductive interface.
            NOVA_REQUIRE(std::abs(center_x - 0.5) > 1.0e-12);
            exact_cell_k.push_back(center_x < 0.5 ? 2.0 : 6.0);
        }

        for (int representation = 0; representation < 2; ++representation) {
            ScalarDiffusion3DProblem problem;
            problem.mesh = box.mesh;
            if (representation == 0) {
                problem.cell_diffusivity = exact_cell_k;
            } else {
                problem.diffusivity = [](const double x, double, double) {
                    return x < 0.5 ? 2.0 : 6.0;
                };
            }

            for (const auto& face : box.boundary_faces) {
                if (all_nodes_on_x(box.mesh, face, 0.0)) {
                    problem.dirichlet_faces.push_back({
                        face, [](double, double, double) { return 0.0; }});
                } else if (all_nodes_on_x(box.mesh, face, 1.0)) {
                    problem.neumann_faces.push_back({
                        face, [](double, double, double) { return -1.0; }});
                }
                // No condition on other faces: physical zero-normal-flux.
            }

            nova::linalg::SolverOptions options;
            options.relative_tolerance = 1.0e-12;
            options.absolute_tolerance = 1.0e-14;
            const auto solved = solve_scalar_diffusion_3d(problem, options);
            NOVA_REQUIRE(solved.linear_solver.converged);
            NOVA_REQUIRE(solved.values.size() == box.mesh.nodes.size());

            for (std::size_t i = 0; i < box.mesh.nodes.size(); ++i) {
                NOVA_REQUIRE_NEAR(
                    solved.values[i],
                    exact_temperature(box.mesh.nodes[i].x), 2.0e-8);
            }

            const auto fields = evaluate_p1_tetrahedral_flux_3d(
                problem, solved.values);
            NOVA_REQUIRE(fields.size() == box.mesh.elements.size());
            double energy = 0.0;
            for (std::size_t i = 0; i < fields.size(); ++i) {
                const double k = exact_cell_k[i];
                NOVA_REQUIRE_NEAR(fields[i].gradient.x, 1.0 / k, 2.0e-7);
                NOVA_REQUIRE_NEAR(fields[i].gradient.y, 0.0, 2.0e-7);
                NOVA_REQUIRE_NEAR(fields[i].gradient.z, 0.0, 2.0e-7);
                NOVA_REQUIRE_NEAR(fields[i].flux.x, -1.0, 4.0e-7);
                NOVA_REQUIRE_NEAR(fields[i].flux.y, 0.0, 4.0e-7);
                NOVA_REQUIRE_NEAR(fields[i].flux.z, 0.0, 4.0e-7);
                const auto& g = fields[i].gradient;
                energy += 0.5 * k * (g.x*g.x + g.y*g.y + g.z*g.z)
                    * quality.elements[i].volume;
            }
            NOVA_REQUIRE_NEAR(energy, 1.0 / 6.0, 2.0e-8);

            // A discontinuous gradient must NOT produce a flux jump when
            // its material-weighted normal component is continuous.
            const auto jumps = evaluate_interior_flux_jumps_3d(
                problem.mesh, fields);
            NOVA_REQUIRE(!jumps.faces.empty());
            NOVA_REQUIRE(jumps.maximum_absolute_normal_jump < 1.0e-6);
            NOVA_REQUIRE(jumps.global_indicator < 1.0e-6);
        }
    }
}
