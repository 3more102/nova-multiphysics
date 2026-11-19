#include "test_framework.hpp"

#include "nova/fem/conservative_face_flux_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using namespace nova::fem;

ScalarDiffusion3DProblem tetra() {
    ScalarDiffusion3DProblem p;
    p.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    p.diffusivity = [](double, double, double) { return 2.0; };
    p.source = [](double, double, double) { return 0.0; };
    return p;
}
bool at_x(const TetrahedralMesh3D& mesh, BoundaryFace3D face, double x) {
    return std::all_of(face.begin(), face.end(), [&](std::size_t id) {
        return mesh.nodes[id].x == x;
    });
}
const ConservativeFaceHeatRate3D& face_at(const ConservativeFaceFlux3DResult& r,
                                         BoundaryFace3D vertices) {
    std::sort(vertices.begin(), vertices.end());
    auto it = std::find_if(r.faces.begin(), r.faces.end(), [&](const auto& row) {
        return row.vertices == vertices;
    });
    NOVA_REQUIRE(it != r.faces.end());
    return *it;
}
template<class F>
bool invalid(F&& f) {
    try { f(); }
    catch (const std::invalid_argument&) { return true; }
    return false;
}
} // namespace

NOVA_TEST(conservative_face_flux_preserves_analytic_affine_heat_rates) {
    const auto p = tetra();
    const std::array u{0.0, 1.0, 0.0, 0.0}; // q = (-2,0,0)
    const auto raw = evaluate_p1_tetrahedral_flux_3d(p, u);
    const auto faces = extract_boundary_faces_3d(p.mesh);
    const auto r = reconstruct_conservative_face_flux_3d(
        p.mesh, raw, std::array{0.0},
        std::span<const PrescribedBoundaryHeatRate3D>{}, faces);
    NOVA_REQUIRE(r.faces.size() == 4);
    NOVA_REQUIRE_NEAR(r.maximum_absolute_cell_defect, 0.0, 1e-12);
    NOVA_REQUIRE_NEAR(face_at(r, {0,2,3}).outward_heat_rate, 1.0, 1e-12);
    NOVA_REQUIRE_NEAR(face_at(r, {1,2,3}).outward_heat_rate, -1.0, 1e-12);
    NOVA_REQUIRE_NEAR(face_at(r, {0,1,2}).outward_heat_rate, 0.0, 1e-12);
}

NOVA_TEST(conservative_face_flux_reconstructs_source_on_fully_fixed_tetra) {
    auto p = tetra();
    p.source = [](double, double, double) { return 60.0; }; // 10W / tetra
    for (std::size_t i = 0; i < 4; ++i) p.dirichlet_constraints.push_back({i, 300.0});
    const auto u = solve_scalar_diffusion_3d(p);
    NOVA_REQUIRE(u.linear_solver.converged);
    const auto r = reconstruct_conservative_face_flux_3d(p, u.values);
    NOVA_REQUIRE(r.faces.size() == 4);
    NOVA_REQUIRE_NEAR(r.maximum_absolute_cell_defect, 0.0, 1e-10);
    double sum = 0.0;
    for (const auto& f : r.faces) sum += f.outward_heat_rate;
    NOVA_REQUIRE_NEAR(sum, 10.0, 1e-10);
}

NOVA_TEST(conservative_face_flux_shared_interior_face_is_equal_and_opposite) {
    auto p = tetra();
    p.mesh.nodes.push_back({0,0,-1});
    p.mesh.elements.push_back({0,2,1,4});
    const std::array<P1TetrahedralFlux3D, 2> raw{{
        {Point3D{}, Point3D{0,0,-1}},
        {Point3D{}, Point3D{0,0,3}}
    }};
    const std::array sources{2.0, -2.0};
    const auto r = reconstruct_conservative_face_flux_3d(
        p.mesh, raw, sources,
        std::span<const PrescribedBoundaryHeatRate3D>{},
        std::span<const BoundaryFace3D>{});
    const auto& shared = face_at(r, {0,1,2});
    NOVA_REQUIRE(shared.second_cell.has_value());
    NOVA_REQUIRE_NEAR(shared.outward_heat_rate, 2.0, 1e-10);
    NOVA_REQUIRE_NEAR(r.maximum_absolute_cell_defect, 0.0, 1e-10);
    for (const auto& f : r.faces) {
        if (!f.second_cell) NOVA_REQUIRE_NEAR(f.outward_heat_rate, 0.0, 1e-12);
    }
}

NOVA_TEST(conservative_face_flux_solves_layered_wall_with_volume_heating_and_robin) {
    auto box = make_box_mesh_3d(0,1,0,1,0,1,2,2,2);
    ScalarDiffusion3DProblem p;
    p.mesh = std::move(box.mesh);
    p.source = [](double, double, double) { return 0.0; };
    for (const auto& tet : p.mesh.elements) {
        double cx = 0;
        for (auto id : tet) cx += p.mesh.nodes[id].x * 0.25;
        p.cell_diffusivity.push_back(cx < 0.5 ? 2.0 : 6.0);
        p.cell_source.push_back(cx < 0.5 ? 120.0 : 360.0);
    }
    for (const auto& f : box.boundary_faces) {
        if (at_x(p.mesh, f, 0)) {
            p.dirichlet_faces.push_back({f, [](double,double,double) {return 400.0;}});
        } else if (at_x(p.mesh, f, 1)) {
            p.robin_faces.push_back({f,
                [](double,double,double) {return 12.0;},
                [](double,double,double) {return 300.0;}});
        }
    }
    const auto solution = solve_scalar_diffusion_3d(p);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto r = reconstruct_conservative_face_flux_3d(p, solution.values);
    NOVA_REQUIRE(r.faces.size() > p.mesh.elements.size());
    NOVA_REQUIRE(r.maximum_absolute_cell_defect < 1e-7);
    double robin=0, dirichlet=0;
    for (const auto& f : r.faces) {
        if (f.second_cell) continue;
        if (at_x(p.mesh, f.vertices, 0)) dirichlet += f.outward_heat_rate;
        else if (at_x(p.mesh, f.vertices, 1)) robin += f.outward_heat_rate;
        else NOVA_REQUIRE_NEAR(f.outward_heat_rate, 0.0, 1e-12);
    }
    NOVA_REQUIRE_NEAR(robin, 384.0, 2e-6);
    NOVA_REQUIRE_NEAR(dirichlet, -144.0, 2e-6);
    NOVA_REQUIRE_NEAR(robin + dirichlet, 240.0, 2e-6);
}

NOVA_TEST(conservative_face_flux_rejects_incompatible_and_duplicate_boundary_roles) {
    const auto p = tetra();
    const std::array<P1TetrahedralFlux3D, 1> raw{{{Point3D{}, Point3D{}}}};
    const auto faces = extract_boundary_faces_3d(p.mesh);
    NOVA_REQUIRE(invalid([&] {
        (void)reconstruct_conservative_face_flux_3d(
            p.mesh, raw, std::array{1.0},
            std::span<const PrescribedBoundaryHeatRate3D>{},
            std::span<const BoundaryFace3D>{});
    }));
    NOVA_REQUIRE(invalid([&] {
        const std::array<BoundaryFace3D, 2> duplicated{faces[0], faces[0]};
        (void)reconstruct_conservative_face_flux_3d(
            p.mesh, raw, std::array{0.0},
            std::span<const PrescribedBoundaryHeatRate3D>{}, duplicated);
    }));
    NOVA_REQUIRE(invalid([&] {
        const std::array<PrescribedBoundaryHeatRate3D, 1> fixed{{{faces[0],1.0}}};
        const std::array<BoundaryFace3D,1> adjustable{faces[0]};
        (void)reconstruct_conservative_face_flux_3d(
            p.mesh, raw, std::array{0.0}, fixed, adjustable);
    }));
    NOVA_REQUIRE(invalid([&] {
        (void)reconstruct_conservative_face_flux_3d(
            p.mesh, raw, std::array{std::numeric_limits<double>::quiet_NaN()},
            std::span<const PrescribedBoundaryHeatRate3D>{}, faces);
    }));
}

NOVA_TEST(conservative_face_flux_checks_closed_components_independently) {
    auto p = tetra();
    p.mesh.nodes.insert(p.mesh.nodes.end(), {
        Point3D{10,0,0}, Point3D{11,0,0},
        Point3D{10,1,0}, Point3D{10,0,1}
    });
    p.mesh.elements.push_back({4,5,6,7});
    const std::array<P1TetrahedralFlux3D,2> raw{{
        {Point3D{}, Point3D{}}, {Point3D{}, Point3D{}}
    }};
    const std::array<double,2> source{1.0e14, 1.0};
    // The giant source in the open component must never mask the 1 W
    // incompatible source in the disconnected, fully insulated component.
    const std::array<BoundaryFace3D,1> open{{{0,1,2}}};
    NOVA_REQUIRE(invalid([&] {
        (void)reconstruct_conservative_face_flux_3d(
            p.mesh, raw, source,
            std::span<const PrescribedBoundaryHeatRate3D>{}, open);
    }));
}
