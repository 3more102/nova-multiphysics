#include "test_framework.hpp"

#include "nova/fem/conservative_heat_flux_3d.hpp"
#include "nova/fem/energy_balance_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using namespace nova::fem;

ScalarDiffusion3DProblem tetra(double source) {
    ScalarDiffusion3DProblem p;
    p.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    p.diffusivity = [](double, double, double) { return 2.0; };
    p.source = [source](double, double, double) { return source; };
    return p;
}

bool on_x(const TetrahedralMesh3D& mesh, const BoundaryFace3D& face, double value) {
    return std::all_of(face.begin(), face.end(), [&](std::size_t node) {
        return mesh.nodes[node].x == value;
    });
}

template<class Action>
bool invalid(Action&& action) {
    try { action(); } catch (const std::invalid_argument&) { return true; }
    return false;
}
} // namespace

NOVA_TEST(conservative_heat_3d_one_tetrahedron_routes_generation_to_dirichlet_face) {
    auto p = tetra(6.0); // 6 W/m^3 * 1/6 m^3 = 1 W
    p.dirichlet_faces.push_back({
        {1,2,3}, [](double, double, double) { return 0.0; }
    });
    const std::array temperature{0.0, 0.0, 0.0, 0.0};
    auto report = reconstruct_conservative_heat_flux_3d(p, temperature);
    NOVA_REQUIRE(report.faces.size() == 4);
    NOVA_REQUIRE_NEAR(report.generated_heat_per_cell[0], 1.0, 1e-13);
    NOVA_REQUIRE_NEAR(report.outward_heat_per_cell[0], 1.0, 1e-13);
    NOVA_REQUIRE(report.max_cell_balance_defect < 1e-12);
    for (const auto& face : report.faces) {
        NOVA_REQUIRE(face.neighbor_cell == ConservativeFaceHeat3D::no_neighbor);
        if (face.vertices == BoundaryFace3D{1,2,3}) {
            NOVA_REQUIRE_NEAR(face.outward_heat_rate, 1.0, 1e-12);
        } else {
            NOVA_REQUIRE_NEAR(face.outward_heat_rate, 0.0, 1e-12);
        }
    }
    std::swap(p.mesh.elements[0][0], p.mesh.elements[0][1]);
    const auto flipped = reconstruct_conservative_heat_flux_3d(p, temperature);
    for (std::size_t i = 0; i < report.faces.size(); ++i) {
        NOVA_REQUIRE(report.faces[i].vertices == flipped.faces[i].vertices);
        NOVA_REQUIRE_NEAR(report.faces[i].outward_heat_rate,
                          flipped.faces[i].outward_heat_rate, 1e-12);
    }
}

NOVA_TEST(conservative_heat_3d_preserves_neumann_flux_and_pure_neumann_gauge) {
    auto p = tetra(6.0);
    p.neumann_faces.push_back({
        {0,2,3}, [](double, double, double) { return 2.0; }
    });
    const std::array temperature{0.0,0.0,0.0,0.0};
    const auto report = reconstruct_conservative_heat_flux_3d(p, temperature);
    NOVA_REQUIRE_NEAR(report.outward_heat_per_cell[0], 1.0, 1e-12);
    NOVA_REQUIRE(report.max_cell_balance_defect < 1e-12);
    for (const auto& f : report.faces) {
        NOVA_REQUIRE_NEAR(f.outward_heat_rate,
                          f.vertices == BoundaryFace3D{0,2,3} ? 1.0 : 0.0,
                          1e-12);
    }
}

NOVA_TEST(conservative_heat_3d_rejects_missing_outlet_or_unbalanced_pure_neumann) {
    auto p = tetra(6.0);
    const std::array temperature{0.0,0.0,0.0,0.0};
    NOVA_REQUIRE(invalid([&] {
        (void)reconstruct_conservative_heat_flux_3d(p, temperature);
    }));
    p.neumann_faces.push_back({
        {0,2,3}, [](double, double, double) { return 1.0; }
    });
    NOVA_REQUIRE(invalid([&] {
        (void)reconstruct_conservative_heat_flux_3d(p, temperature);
    }));
    p.neumann_faces.clear();
    p.dirichlet_constraints.push_back({0,0.0});
    NOVA_REQUIRE(invalid([&] {
        (void)reconstruct_conservative_heat_flux_3d(p, temperature);
    }));
}

NOVA_TEST(conservative_heat_3d_anisotropic_affine_patch_preserves_exact_face_rates) {
    auto p = tetra(0.0);
    p.diffusivity = {};
    p.diffusivity_tensor = [](double, double, double) {
        return SymmetricTensor3D{2.0,3.0,4.0,0.0,0.0,0.0};
    };
    for (const auto face : extract_boundary_faces_3d(p.mesh)) {
        p.dirichlet_faces.push_back({
            face, [](double x, double, double) { return x; }
        });
    }
    const std::array temperature{0.0,1.0,0.0,0.0};
    const auto report = reconstruct_conservative_heat_flux_3d(p, temperature);
    NOVA_REQUIRE(report.max_cell_balance_defect < 1e-12);
    for (const auto& f : report.faces) {
        if (f.vertices == BoundaryFace3D{0,2,3}) {
            NOVA_REQUIRE_NEAR(f.outward_heat_rate, 1.0, 1e-12);
        } else if (f.vertices == BoundaryFace3D{1,2,3}) {
            NOVA_REQUIRE_NEAR(f.outward_heat_rate, -1.0, 1e-12);
        } else {
            NOVA_REQUIRE_NEAR(f.outward_heat_rate, 0.0, 1e-12);
        }
    }
}

NOVA_TEST(conservative_heat_3d_multicell_anisotropic_source_and_robin_balance) {
    auto box = make_box_mesh_3d(0,1,0,1,0,1,2,2,2);
    ScalarDiffusion3DProblem p;
    p.mesh = std::move(box.mesh);
    p.source = [](double, double, double) { return 0.0; };
    for (const auto& tetra : p.mesh.elements) {
        double xc = 0;
        for (const auto node : tetra) xc += p.mesh.nodes[node].x / 4.0;
        p.cell_diffusivity_tensors.push_back({
            xc < 0.5 ? 2.0 : 6.0, 4.0, 5.0, 0.0, 0.0, 0.0
        });
        p.cell_source.push_back(xc < 0.5 ? 120.0 : 360.0);
    }
    for (const auto& face : box.boundary_faces) {
        if (on_x(p.mesh, face, 0.0)) {
            p.dirichlet_faces.push_back({face, [](double, double, double) {
                return 400.0;
            }});
        } else if (on_x(p.mesh, face, 1.0)) {
            p.robin_faces.push_back({
                face,
                [](double, double, double) { return 12.0; },
                [](double, double, double) { return 300.0; }
            });
        }
    }
    const auto solved = solve_scalar_diffusion_3d(p);
    NOVA_REQUIRE(solved.linear_solver.converged);
    const auto audit = evaluate_diffusion_energy_balance_3d(p, solved.values);
    const auto report = reconstruct_conservative_heat_flux_3d(p, solved.values);
    NOVA_REQUIRE(report.generated_heat_per_cell.size() == p.mesh.elements.size());
    NOVA_REQUIRE_NEAR(audit.generated_heat, 240.0, 2e-9);
    NOVA_REQUIRE(report.max_cell_balance_defect < 2e-7);
    std::size_t shared = 0;
    double total_exterior = 0;
    double robin = 0;
    for (const auto& f : report.faces) {
        if (f.neighbor_cell != ConservativeFaceHeat3D::no_neighbor) {
            ++shared;
            NOVA_REQUIRE(f.owner_cell != f.neighbor_cell);
            continue;
        }
        total_exterior += f.outward_heat_rate;
        if (on_x(p.mesh, f.vertices, 1.0)) robin += f.outward_heat_rate;
    }
    NOVA_REQUIRE(shared > 0);
    NOVA_REQUIRE_NEAR(total_exterior, 240.0, 2e-7);
    NOVA_REQUIRE_NEAR(robin, audit.robin_outward_heat, 2e-8);
    for (std::size_t i = 0; i < p.mesh.elements.size(); ++i) {
        NOVA_REQUIRE_NEAR(report.outward_heat_per_cell[i],
                          report.generated_heat_per_cell[i], 2e-7);
    }
}
