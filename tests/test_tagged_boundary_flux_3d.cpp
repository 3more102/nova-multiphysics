#include "test_framework.hpp"

#include "nova/fem/tagged_boundary_flux_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using namespace nova::fem;

ScalarDiffusion3DProblem unit_problem() {
    ScalarDiffusion3DProblem problem;
    problem.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    problem.diffusivity = [](double, double, double) { return 2.0; };
    return problem;
}

int classify(const TetrahedralMesh3D& mesh, const BoundaryFace3D& face) {
    bool x_zero = true;
    bool far_face = true;
    for (std::size_t node : face) {
        const auto& p = mesh.nodes[node];
        x_zero = x_zero && p.x == 0.0;
        far_face = far_face && (p.x + p.y + p.z) == 1.0;
    }
    if (x_zero) return 101;
    if (far_face) return 102;
    return 0;
}

const TaggedBoundaryFlux3D& group(const std::vector<TaggedBoundaryFlux3D>& rows, int tag) {
    const auto found = std::find_if(rows.begin(), rows.end(),
                                    [tag](const auto& r) { return r.physical_tag == tag; });
    NOVA_REQUIRE(found != rows.end());
    return *found;
}

template<class F>
bool rejects(F&& operation) {
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

} // namespace

NOVA_TEST(tagged_boundary_flux_reports_signed_analytic_heat_rates) {
    auto problem = unit_problem();
    const auto faces = extract_boundary_faces_3d(problem.mesh);
    std::vector<int> tags;
    for (const auto& face : faces) tags.push_back(classify(problem.mesh, face));
    const std::array nodal{0.0, 1.0, 0.0, 0.0}; // T = x, q = (-2,0,0)
    const auto field = evaluate_p1_tetrahedral_flux_3d(problem, nodal);
    const auto result = evaluate_tagged_boundary_flux_3d(
        problem.mesh, faces, tags, field
    );

    NOVA_REQUIRE(result.size() == 3);
    NOVA_REQUIRE(result[0].physical_tag == 0);
    NOVA_REQUIRE(group(result, 0).face_count == 2);
    NOVA_REQUIRE_NEAR(group(result, 0).area, 1.0, 1e-14);
    NOVA_REQUIRE_NEAR(group(result, 0).outward_heat_rate, 0.0, 1e-14);
    NOVA_REQUIRE(group(result, 101).face_count == 1);
    NOVA_REQUIRE_NEAR(group(result, 101).area, 0.5, 1e-14);
    NOVA_REQUIRE_NEAR(group(result, 101).outward_heat_rate, 1.0, 1e-14);
    NOVA_REQUIRE(group(result, 102).face_count == 1);
    NOVA_REQUIRE_NEAR(group(result, 102).area, std::sqrt(3.0) / 2, 1e-14);
    NOVA_REQUIRE_NEAR(group(result, 102).outward_heat_rate, -1.0, 1e-14);
}

NOVA_TEST(tagged_boundary_flux_ignores_boundary_winding_and_cell_orientation) {
    auto problem = unit_problem();
    auto faces = extract_boundary_faces_3d(problem.mesh);
    std::vector<int> tags;
    for (const auto& face : faces) tags.push_back(classify(problem.mesh, face));
    std::reverse(faces.begin(), faces.end());
    std::reverse(tags.begin(), tags.end());
    for (auto& face : faces) std::swap(face[0], face[2]);
    std::swap(problem.mesh.elements[0][0], problem.mesh.elements[0][1]);
    const auto fields = evaluate_p1_tetrahedral_flux_3d(
        problem, std::array{0.0, 1.0, 0.0, 0.0}
    );
    const auto rows = evaluate_tagged_boundary_flux_3d(
        problem.mesh, faces, tags, fields
    );
    NOVA_REQUIRE_NEAR(group(rows, 101).outward_heat_rate, 1.0, 1e-14);
    NOVA_REQUIRE_NEAR(group(rows, 102).outward_heat_rate, -1.0, 1e-14);
}

NOVA_TEST(tagged_boundary_flux_rejects_duplicate_missing_interior_and_invalid_data) {
    auto problem = unit_problem();
    const auto faces = extract_boundary_faces_3d(problem.mesh);
    const auto fields = evaluate_p1_tetrahedral_flux_3d(
        problem, std::array{0.0, 1.0, 0.0, 0.0}
    );
    auto valid = [&] {
        return evaluate_tagged_boundary_flux_3d(
            problem.mesh, faces, std::array{0, 0, 0, 0}, fields);
    };
    NOVA_REQUIRE(valid().size() == 1);
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_tagged_boundary_flux_3d(
            problem.mesh, faces, std::array{0, 0, 0}, fields);
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_tagged_boundary_flux_3d(
            problem.mesh, std::array{faces[0], faces[0]}, std::array{7, 7}, fields);
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_tagged_boundary_flux_3d(
            problem.mesh, std::array<BoundaryFace3D, 1>{faces[0]}, std::array{-1}, fields);
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_tagged_boundary_flux_3d(
            problem.mesh, std::array<BoundaryFace3D, 1>{BoundaryFace3D{0,1,99}}, std::array{1}, fields);
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_tagged_boundary_flux_3d(
            problem.mesh, faces, std::array{0,0,0,0},
            std::array<P1TetrahedralFlux3D, 0>{});
    }));
    auto nan_field = fields;
    nan_field[0].flux.y = std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_tagged_boundary_flux_3d(
            problem.mesh, faces, std::array{0,0,0,0}, nan_field);
    }));

    // The triangular face z=0 is shared by two tetrahedra, and is not
    // available for a physical exterior heat-flow measurement.
    problem.mesh.nodes.push_back({0,0,-1});
    problem.mesh.elements.push_back({0,1,2,4});
    const std::array two_fields{
        P1TetrahedralFlux3D{Point3D{}, Point3D{-2,0,0}},
        P1TetrahedralFlux3D{Point3D{}, Point3D{-2,0,0}}
    };
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_tagged_boundary_flux_3d(
            problem.mesh, std::array<BoundaryFace3D, 1>{BoundaryFace3D{0,1,2}}, std::array{9}, two_fields);
    }));
}
