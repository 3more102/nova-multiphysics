#include "test_framework.hpp"

#include "nova/fem/interior_flux_jumps_3d.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using namespace nova::fem;

TetrahedralMesh3D adjacent_tetrahedra() {
    // Shared z=0 face: area=1/2 and diameter=sqrt(2).
    return {
        {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {0,0,-1}},
        {{0,1,2,3}, {0,2,1,4}}
    };
}

template<class Function>
bool rejects_invalid(Function&& action) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}
} // namespace

NOVA_TEST(jumps_3d_zero_for_constant_flux_and_one_tetrahedron) {
    const auto mesh = adjacent_tetrahedra();
    const std::array<P1TetrahedralFlux3D, 2> fields{{
        {{1,2,3}, {2,-4,6}},
        {{5,7,8}, {2,-4,6}}
    }};
    const auto result = evaluate_interior_flux_jumps_3d(mesh, fields);
    NOVA_REQUIRE(result.faces.size() == 1U);
    NOVA_REQUIRE(result.cell_indicators.size() == 2U);
    NOVA_REQUIRE_NEAR(result.faces[0].area, 0.5, 1e-15);
    NOVA_REQUIRE_NEAR(result.faces[0].diameter, std::sqrt(2.0), 1e-15);
    NOVA_REQUIRE_NEAR(result.faces[0].normal_jump, 0, 1e-15);
    NOVA_REQUIRE_NEAR(result.global_indicator, 0, 1e-15);
    NOVA_REQUIRE_NEAR(result.maximum_absolute_normal_jump, 0, 1e-15);
    NOVA_REQUIRE_NEAR(result.cell_indicators[0], 0, 1e-15);
    NOVA_REQUIRE_NEAR(result.cell_indicators[1], 0, 1e-15);

    TetrahedralMesh3D isolated = mesh;
    isolated.elements.pop_back();
    isolated.nodes.pop_back();
    const std::array<P1TetrahedralFlux3D, 1> one{{fields[0]}};
    const auto alone = evaluate_interior_flux_jumps_3d(isolated, one);
    NOVA_REQUIRE(alone.faces.empty());
    NOVA_REQUIRE(alone.cell_indicators.size() == 1U);
    NOVA_REQUIRE_NEAR(alone.global_indicator, 0, 0);
}

NOVA_TEST(jumps_3d_exact_face_weight_and_cell_partition) {
    const auto mesh = adjacent_tetrahedra();
    const std::array<P1TetrahedralFlux3D, 2> fields{{
        {{}, {0,0,2}}, {{}, {0,0,4}}
    }};
    const auto result = evaluate_interior_flux_jumps_3d(mesh, fields);
    const double contribution = 2.0 * std::sqrt(2.0); // h_F * |F| * (2)^2
    NOVA_REQUIRE_NEAR(result.faces[0].normal_jump, 2.0, 1e-14);
    NOVA_REQUIRE_NEAR(result.maximum_absolute_normal_jump, 2.0, 1e-14);
    NOVA_REQUIRE_NEAR(result.global_indicator, std::sqrt(contribution), 1e-14);
    for (const double cell : result.cell_indicators) {
        NOVA_REQUIRE_NEAR(cell, std::sqrt(0.5 * contribution), 1e-14);
    }
    NOVA_REQUIRE_NEAR(
        result.cell_indicators[0]*result.cell_indicators[0] +
        result.cell_indicators[1]*result.cell_indicators[1],
        result.global_indicator*result.global_indicator, 1e-14
    );
}

NOVA_TEST(jumps_3d_invariant_to_cell_and_vertex_reordering) {
    auto mesh = adjacent_tetrahedra();
    std::array<P1TetrahedralFlux3D, 2> fields{{
        {{}, {0,0,2}}, {{}, {0,0,4}}
    }};
    const auto reference = evaluate_interior_flux_jumps_3d(mesh, fields);
    std::swap(mesh.elements[0], mesh.elements[1]);
    std::swap(fields[0], fields[1]);
    std::swap(mesh.elements[0][0], mesh.elements[0][1]);
    std::swap(mesh.elements[1][1], mesh.elements[1][2]);
    const auto reordered = evaluate_interior_flux_jumps_3d(mesh, fields);
    NOVA_REQUIRE_NEAR(reordered.global_indicator, reference.global_indicator, 1e-14);
    NOVA_REQUIRE_NEAR(reordered.faces[0].normal_jump,
                      reference.faces[0].normal_jump, 1e-14);
}

NOVA_TEST(jumps_3d_detect_material_interface_normal_flux_mismatch) {
    ScalarDiffusion3DProblem problem;
    problem.mesh = adjacent_tetrahedra();
    problem.cell_diffusivity = {2.0, 6.0};
    // Piecewise affine temperature is continuous across z=0, with
    // k_top * dT_top/dz = k_bottom * dT_bottom/dz = 6.
    const std::array compatible{0.0, 0.0, 0.0, 3.0, -1.0};
    const auto exact_flux = evaluate_p1_tetrahedral_flux_3d(problem, compatible);
    const auto matched = evaluate_interior_flux_jumps_3d(problem.mesh, exact_flux);
    NOVA_REQUIRE_NEAR(matched.global_indicator, 0, 1e-13);

    // Same slope on both sides yields a nonphysical flux jump of magnitude 4.
    const std::array mismatched{0.0, 0.0, 0.0, 1.0, -1.0};
    const auto bad_flux = evaluate_p1_tetrahedral_flux_3d(problem, mismatched);
    const auto bad = evaluate_interior_flux_jumps_3d(problem.mesh, bad_flux);
    NOVA_REQUIRE_NEAR(bad.maximum_absolute_normal_jump, 4.0, 1e-13);
    NOVA_REQUIRE(bad.global_indicator > 0.0);
}

NOVA_TEST(jumps_3d_reject_nonfinite_or_wrong_sized_fluxes) {
    const auto mesh = adjacent_tetrahedra();
    const std::array<P1TetrahedralFlux3D, 1> short_fields{{{{},{0,0,1}}}};
    NOVA_REQUIRE(rejects_invalid([&] {
        (void)evaluate_interior_flux_jumps_3d(mesh, short_fields);
    }));
    std::array<P1TetrahedralFlux3D, 2> fields{{{{},{0,0,1}},{{},{0,0,2}}}};
    fields[1].flux.z = std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects_invalid([&] {
        (void)evaluate_interior_flux_jumps_3d(mesh, fields);
    }));
}
