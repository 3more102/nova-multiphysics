#include "test_framework.hpp"

#include "nova/fem/conservative_face_heat_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using namespace nova::fem;

TetrahedralMesh3D two_tetrahedra() {
    return {
        {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {0,0,-1}},
        {{0,1,2,3}, {0,2,1,4}}
    };
}

ScalarDiffusion3DProblem example_problem() {
    ScalarDiffusion3DProblem p;
    p.mesh = two_tetrahedra();
    p.cell_diffusivity = {2.0, 6.0};
    p.source = [](double, double, double) { return 0.0; };
    return p;
}

template<class F>
bool rejects_invalid(F&& callback) {
    try {
        callback();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

void fix_all_exterior(ScalarDiffusion3DProblem& p, double value) {
    for (const auto& face : extract_boundary_faces_3d(p.mesh)) {
        p.dirichlet_faces.push_back({
            face, [value](double, double, double) { return value; }
        });
    }
}
} // namespace

NOVA_TEST(face_heat_3d_exact_two_material_flux_patch_has_no_correction) {
    auto p = example_problem();
    for (const auto& face : extract_boundary_faces_3d(p.mesh)) {
        p.dirichlet_faces.push_back({
            face, [](double, double, double z) { return z >= 0.0 ? 3.0*z : z; }
        });
    }
    const std::array values{0.0, 0.0, 0.0, 3.0, -1.0};
    const auto report = reconcile_p1_face_heat_rates_3d(p, values);
    NOVA_REQUIRE(report.faces.size() == 7U);
    NOVA_REQUIRE(report.floating_components == 0U);
    NOVA_REQUIRE(report.multiplier_solver.converged);
    NOVA_REQUIRE_NEAR(report.total_generated_heat, 0.0, 1e-14);
    NOVA_REQUIRE_NEAR(report.total_exterior_outward_heat, 0.0, 1e-12);
    NOVA_REQUIRE_NEAR(report.max_absolute_correction, 0.0, 1e-11);
    NOVA_REQUIRE_NEAR(report.max_absolute_cell_defect, 0.0, 1e-11);
    std::size_t count_internal = 0;
    for (const auto& f : report.faces) {
        if (f.interior) {
            ++count_internal;
            NOVA_REQUIRE(f.condition == FaceHeatCondition3D::Interior);
            NOVA_REQUIRE_NEAR(f.area, 0.5, 1e-14);
            NOVA_REQUIRE_NEAR(f.outward_heat, 3.0, 1e-12);
        }
    }
    NOVA_REQUIRE(count_internal == 1U);
}

NOVA_TEST(face_heat_3d_recovers_nonzero_cell_generation_and_local_balances) {
    auto p = example_problem();
    p.source = [](double, double, double) { return 6.0; };
    fix_all_exterior(p, 0.0);
    const std::array values{0.0, 0.0, 0.0, 0.0, 0.0};
    const auto report = reconcile_p1_face_heat_rates_3d(p, values);
    NOVA_REQUIRE_NEAR(report.total_generated_heat, 2.0, 2e-13);
    NOVA_REQUIRE_NEAR(report.total_exterior_outward_heat, 2.0, 2e-10);
    NOVA_REQUIRE_NEAR(report.generated_heat_per_cell[0], 1.0, 1e-13);
    NOVA_REQUIRE_NEAR(report.generated_heat_per_cell[1], 1.0, 1e-13);
    NOVA_REQUIRE(report.max_absolute_correction > 0.0);
    NOVA_REQUIRE(report.max_absolute_cell_defect < 1e-9);
    for (double defect : report.cell_balance_defect) {
        NOVA_REQUIRE_NEAR(defect, 0.0, 1e-9);
    }
}

NOVA_TEST(face_heat_3d_preserves_neumann_and_robin_heat_exactly) {
    auto p = example_problem();
    const auto neumann = BoundaryFace3D{0,1,3};
    const auto robin = BoundaryFace3D{0,1,4};
    for (auto face : extract_boundary_faces_3d(p.mesh)) {
        auto sorted = face;
        std::sort(sorted.begin(), sorted.end());
        if (sorted == neumann) {
            p.neumann_faces.push_back({
                face, [](double, double, double) { return 2.0; }
            });
        } else if (sorted == robin) {
            p.robin_faces.push_back({
                face,
                [](double, double, double) { return 2.0; },
                [](double, double, double) { return 5.0; }
            });
        } else {
            p.dirichlet_faces.push_back({
                face, [](double, double, double) { return 0.0; }
            });
        }
    }
    const std::array values{0.0,0.0,0.0,0.0,0.0};
    const auto report = reconcile_p1_face_heat_rates_3d(p, values);
    bool found_neumann = false;
    bool found_robin = false;
    for (const auto& face : report.faces) {
        if (face.condition == FaceHeatCondition3D::Neumann) {
            found_neumann = true;
            NOVA_REQUIRE_NEAR(face.outward_heat, 1.0, 1e-13);
        } else if (face.condition == FaceHeatCondition3D::Robin) {
            found_robin = true;
            NOVA_REQUIRE_NEAR(face.outward_heat, -5.0, 1e-13);
        }
    }
    NOVA_REQUIRE(found_neumann && found_robin);
    NOVA_REQUIRE_NEAR(report.total_exterior_outward_heat, 0.0, 1e-9);
    NOVA_REQUIRE(report.max_absolute_cell_defect < 1e-9);
}

NOVA_TEST(face_heat_3d_orientation_and_winding_independent) {
    auto p = example_problem();
    p.source = [](double, double, double) { return 6.0; };
    fix_all_exterior(p, 0.0);
    const std::array values{0.0,0.0,0.0,0.0,0.0};
    const auto before = reconcile_p1_face_heat_rates_3d(p, values);
    std::swap(p.mesh.elements[0][0], p.mesh.elements[0][1]);
    std::swap(p.mesh.elements[1][1], p.mesh.elements[1][2]);
    const auto after = reconcile_p1_face_heat_rates_3d(p, values);
    NOVA_REQUIRE(after.faces.size() == before.faces.size());
    for (std::size_t i = 0; i < before.faces.size(); ++i) {
        NOVA_REQUIRE(after.faces[i].nodes == before.faces[i].nodes);
        NOVA_REQUIRE_NEAR(after.faces[i].outward_heat,
                          before.faces[i].outward_heat, 1e-11);
    }
}

NOVA_TEST(face_heat_3d_heterogeneous_spd_cell_tensors_and_sources) {
    auto p = example_problem();
    const auto box = make_box_mesh_3d(0,2,0,1,0,1,2,1,1);
    p.mesh = box.mesh;
    p.cell_diffusivity.clear();
    p.cell_diffusivity_tensors.clear();
    p.cell_source.clear();
    for (std::size_t i = 0; i < p.mesh.elements.size(); ++i) {
        p.cell_diffusivity_tensors.push_back(
            i % 2 == 0 ? SymmetricTensor3D{2,3,4,0.2,0,0} :
                         SymmetricTensor3D{5,3,2,-0.2,0.1,0}
        );
        p.cell_source.push_back(i % 2 == 0 ? 10.0 : 20.0);
    }
    p.source = [](double, double, double) { return 1.0; };
    for (const auto& face : box.boundary_faces) {
        bool fixed = true;
        const double x = p.mesh.nodes[face[0]].x;
        for (const auto node : face) {
            fixed = fixed && p.mesh.nodes[node].x == x;
        }
        if (fixed && (x == 0.0 || x == 2.0)) {
            p.dirichlet_faces.push_back({
                face, [](double, double, double) { return 0.0; }
            });
        }
    }
    std::vector<double> values(p.mesh.nodes.size(), 0.0);
    const auto report = reconcile_p1_face_heat_rates_3d(p, values);
    NOVA_REQUIRE_NEAR(report.total_generated_heat, 32.0, 1e-11);
    NOVA_REQUIRE_NEAR(report.total_exterior_outward_heat, 32.0, 1e-7);
    NOVA_REQUIRE(report.max_absolute_cell_defect < 1e-7);
    NOVA_REQUIRE(report.multiplier_solver.converged);
}

NOVA_TEST(face_heat_3d_pure_neumann_components_require_compatible_sources) {
    auto p = example_problem();
    const std::array values{0.0,0.0,0.0,0.0,0.0};
    const auto balanced = reconcile_p1_face_heat_rates_3d(p, values);
    NOVA_REQUIRE(balanced.floating_components == 1U);
    NOVA_REQUIRE_NEAR(balanced.total_exterior_outward_heat, 0.0, 0);
    NOVA_REQUIRE_NEAR(balanced.max_absolute_cell_defect, 0.0, 1e-13);
    p.source = [](double, double, double) { return 6.0; };
    NOVA_REQUIRE(rejects_invalid([&] {
        (void)reconcile_p1_face_heat_rates_3d(p, values);
    }));
}

NOVA_TEST(face_heat_3d_rejects_nodal_reactions_and_invalid_inputs) {
    auto p = example_problem();
    fix_all_exterior(p, 0.0);
    const std::array values{0.0,0.0,0.0,0.0,0.0};
    p.dirichlet_constraints.push_back({0,0.0});
    NOVA_REQUIRE(rejects_invalid([&] {
        (void)reconcile_p1_face_heat_rates_3d(p, values);
    }));
    p.dirichlet_constraints.clear();
    NOVA_REQUIRE(rejects_invalid([&] {
        (void)reconcile_p1_face_heat_rates_3d(p, std::array{0.0, 1.0});
    }));
    const std::array bad{
        0.0,0.0,0.0,0.0,std::numeric_limits<double>::quiet_NaN()
    };
    NOVA_REQUIRE(rejects_invalid([&] {
        (void)reconcile_p1_face_heat_rates_3d(p, bad);
    }));
    p.dirichlet_faces.push_back({
        {0,1,2}, [](double, double, double) { return 0.0; }
    });
    NOVA_REQUIRE(rejects_invalid([&] {
        (void)reconcile_p1_face_heat_rates_3d(p, values);
    }));
}
