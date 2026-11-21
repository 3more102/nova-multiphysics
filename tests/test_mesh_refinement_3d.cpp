#include "test_framework.hpp"

#include "nova/fem/mesh_refinement_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using nova::fem::BoundaryFace3D;
using nova::fem::TetrahedralMesh3D;
using FaceKey = std::array<std::size_t, 3>;

FaceKey canonical(FaceKey face) {
    std::sort(face.begin(), face.end());
    return face;
}

bool same_outward_winding(const BoundaryFace3D& left,
                          const BoundaryFace3D& right) {
    return left == right ||
        (left[0] == right[1] && left[1] == right[2] && left[2] == right[0]) ||
        (left[0] == right[2] && left[1] == right[0] && left[2] == right[1]);
}

template <class F>
bool throws_mesh_error(F&& action) {
    try {
        action();
    } catch (const std::exception&) {
        return true;
    }
    return false;
}

std::map<FaceKey, std::size_t> face_incidence(const TetrahedralMesh3D& mesh) {
    std::map<FaceKey, std::size_t> counts;
    for (const auto& tet : mesh.elements) {
        for (const auto face : std::array<FaceKey, 4>{{
            {{tet[0], tet[1], tet[2]}},
            {{tet[0], tet[1], tet[3]}},
            {{tet[0], tet[2], tet[3]}},
            {{tet[1], tet[2], tet[3]}}
        }}) {
            ++counts[canonical(face)];
        }
    }
    return counts;
}

} // namespace

NOVA_TEST(tetra_uniform_refinement_is_conforming_and_preserves_volume) {
    const auto coarse = nova::fem::make_box_mesh_3d(
        -1.0, 2.0, 0.0, 2.0, -0.5, 1.5, 2U, 2U, 1U
    );
    const auto fine = nova::fem::refine_tetrahedra_uniformly_3d(
        coarse.mesh, coarse.boundary_faces
    );
    fine.mesh.validate();

    std::set<std::array<std::size_t, 2>> edges;
    for (const auto& tet : coarse.mesh.elements) {
        for (std::size_t a = 0; a < 4U; ++a) {
            for (std::size_t b = a + 1U; b < 4U; ++b) {
                edges.insert({std::min(tet[a], tet[b]), std::max(tet[a], tet[b])});
            }
        }
    }

    NOVA_REQUIRE(fine.coarse_node_count == coarse.mesh.nodes.size());
    NOVA_REQUIRE(fine.coarse_element_count == coarse.mesh.elements.size());
    NOVA_REQUIRE(fine.mesh.nodes.size() == coarse.mesh.nodes.size() + edges.size());
    NOVA_REQUIRE(fine.midpoint_parent_edges.size() == edges.size());
    NOVA_REQUIRE(fine.mesh.elements.size() == 8U * coarse.mesh.elements.size());
    NOVA_REQUIRE(fine.parent_elements.size() == fine.mesh.elements.size());
    NOVA_REQUIRE(fine.boundary_faces.size() == 4U * coarse.boundary_faces.size());
    NOVA_REQUIRE(fine.parent_boundary_faces.size() == fine.boundary_faces.size());

    const auto quality0 = nova::fem::evaluate_mesh_quality_3d(coarse.mesh);
    const auto quality1 = nova::fem::evaluate_mesh_quality_3d(fine.mesh);
    NOVA_REQUIRE_NEAR(quality0.total_volume, quality1.total_volume, 2.0e-12);
    NOVA_REQUIRE(quality1.minimum_mean_ratio > 0.0);
    for (std::size_t child = 0; child < fine.mesh.elements.size(); ++child) {
        NOVA_REQUIRE(fine.parent_elements[child] == child / 8U);
    }

    const auto expected = nova::fem::extract_boundary_faces_3d(fine.mesh);
    NOVA_REQUIRE(expected.size() == fine.boundary_faces.size());
    std::map<FaceKey, BoundaryFace3D> expected_by_key;
    for (const auto face : expected) {
        NOVA_REQUIRE(expected_by_key.emplace(canonical(face), face).second);
    }
    for (std::size_t i = 0; i < fine.boundary_faces.size(); ++i) {
        const auto face = fine.boundary_faces[i];
        const auto found = expected_by_key.find(canonical(face));
        NOVA_REQUIRE(found != expected_by_key.end());
        NOVA_REQUIRE(same_outward_winding(face, found->second));
        NOVA_REQUIRE(fine.parent_boundary_faces[i] == i / 4U);
    }
    const auto incidence = face_incidence(fine.mesh);
    for (const auto& [face, count] : incidence) {
        (void)face;
        NOVA_REQUIRE(count == 1U || count == 2U);
    }
}

NOVA_TEST(tetra_uniform_refinement_transfers_affine_p1_and_cellwise_data) {
    const auto coarse = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 2U, 1U, 1U
    );
    const auto fine = nova::fem::refine_tetrahedra_uniformly_3d(coarse.mesh);
    std::vector<double> coarse_temperature;
    for (const auto& node : coarse.mesh.nodes) {
        coarse_temperature.push_back(5.0 + 2.0 * node.x
            - 3.0 * node.y + 4.0 * node.z);
    }
    const auto fine_temperature = nova::fem::prolongate_p1_uniform_3d(
        fine, coarse_temperature
    );
    NOVA_REQUIRE(fine_temperature.size() == fine.mesh.nodes.size());
    for (std::size_t i = 0; i < fine.mesh.nodes.size(); ++i) {
        const auto& node = fine.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(fine_temperature[i],
            5.0 + 2.0 * node.x - 3.0 * node.y + 4.0 * node.z, 2.0e-13);
    }

    std::vector<double> sources(coarse.mesh.elements.size());
    for (std::size_t i = 0; i < sources.size(); ++i) {
        sources[i] = 10.0 + static_cast<double>(i);
    }
    const auto child_sources =
        nova::fem::prolongate_cellwise_constant_uniform_3d(fine, sources);
    const auto q0 = nova::fem::evaluate_mesh_quality_3d(coarse.mesh);
    const auto q1 = nova::fem::evaluate_mesh_quality_3d(fine.mesh);
    double coarse_integral = 0.0;
    double refined_integral = 0.0;
    for (std::size_t i = 0; i < sources.size(); ++i) {
        coarse_integral += sources[i] * q0.elements[i].volume;
    }
    for (std::size_t i = 0; i < child_sources.size(); ++i) {
        NOVA_REQUIRE(child_sources[i] == sources[fine.parent_elements[i]]);
        refined_integral += child_sources[i] * q1.elements[i].volume;
    }
    NOVA_REQUIRE_NEAR(refined_integral, coarse_integral, 1.0e-12);
    NOVA_REQUIRE(fine.boundary_faces.empty());

    // Repeated global refinement must preserve the same affine field.
    const auto finer = nova::fem::refine_tetrahedra_uniformly_3d(fine.mesh);
    const auto twice = nova::fem::prolongate_p1_uniform_3d(
        finer, fine_temperature
    );
    for (std::size_t i = 0; i < twice.size(); ++i) {
        const auto& node = finer.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(twice[i],
            5.0 + 2.0 * node.x - 3.0 * node.y + 4.0 * node.z, 3.0e-13);
    }
}

NOVA_TEST(tetra_uniform_refinement_accepts_reversed_parent_orientation) {
    TetrahedralMesh3D tetra{
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
         {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}},
        {{0U, 2U, 1U, 3U}}
    };
    tetra.validate();
    auto boundary = nova::fem::extract_boundary_faces_3d(tetra);
    for (auto& face : boundary) {
        std::rotate(face.begin(), face.begin() + 1, face.end());
    }
    const auto fine = nova::fem::refine_tetrahedra_uniformly_3d(tetra, boundary);
    NOVA_REQUIRE(fine.mesh.elements.size() == 8U);
    NOVA_REQUIRE(fine.boundary_faces.size() == 16U);
    NOVA_REQUIRE_NEAR(
        nova::fem::evaluate_mesh_quality_3d(fine.mesh).total_volume,
        1.0 / 6.0, 1.0e-14
    );
}

NOVA_TEST(tetra_uniform_refinement_rejects_invalid_boundary_faces) {
    const auto coarse = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U
    );
    const auto face = coarse.boundary_faces.front();
    NOVA_REQUIRE(throws_mesh_error([&] {
        nova::fem::refine_tetrahedra_uniformly_3d(
            coarse.mesh, {face, face}
        );
    }));
    NOVA_REQUIRE(throws_mesh_error([&] {
        nova::fem::refine_tetrahedra_uniformly_3d(
            coarse.mesh, {{face[0], face[2], face[1]}}
        );
    }));
    NOVA_REQUIRE(throws_mesh_error([&] {
        nova::fem::refine_tetrahedra_uniformly_3d(
            coarse.mesh, {{coarse.mesh.nodes.size(), face[1], face[2]}}
        );
    }));
    const auto incidence = face_incidence(coarse.mesh);
    bool checked_interior = false;
    for (const auto& [interior, count] : incidence) {
        if (count == 2U) {
            checked_interior = true;
            NOVA_REQUIRE(throws_mesh_error([&] {
                nova::fem::refine_tetrahedra_uniformly_3d(
                    coarse.mesh, {interior}
                );
            }));
            break;
        }
    }
    NOVA_REQUIRE(checked_interior);
}

NOVA_TEST(tetra_uniform_refinement_rejects_invalid_transfer_vectors) {
    const auto coarse = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U
    );
    const auto fine = nova::fem::refine_tetrahedra_uniformly_3d(coarse.mesh);
    std::vector<double> field(coarse.mesh.nodes.size(), 1.0);
    NOVA_REQUIRE(throws_mesh_error([&] {
        nova::fem::prolongate_p1_uniform_3d(fine, std::vector<double>{1.0});
    }));
    field[0] = std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(throws_mesh_error([&] {
        nova::fem::prolongate_p1_uniform_3d(fine, field);
    }));
    std::vector<double> cells(coarse.mesh.elements.size(), 1.0);
    NOVA_REQUIRE(throws_mesh_error([&] {
        nova::fem::prolongate_cellwise_constant_uniform_3d(
            fine, std::vector<double>{1.0}
        );
    }));
    cells[0] = std::numeric_limits<double>::infinity();
    NOVA_REQUIRE(throws_mesh_error([&] {
        nova::fem::prolongate_cellwise_constant_uniform_3d(fine, cells);
    }));
    cells[0] = 1.0;
    auto corrupted = fine;
    corrupted.parent_elements[0] = coarse.mesh.elements.size();
    NOVA_REQUIRE(throws_mesh_error([&] {
        nova::fem::prolongate_cellwise_constant_uniform_3d(corrupted, cells);
    }));
}


NOVA_TEST(tetra_dorfler_marking_selects_the_minimal_deterministic_bulk_set) {
    const std::vector<double> indicators{1.0, 2.0, 3.0, 4.0};
    const auto half = nova::fem::mark_dorfler_elements_3d(indicators, 0.5);
    NOVA_REQUIRE(half.size() == 1U);
    NOVA_REQUIRE(half[0] == 3U);

    const auto eighty_percent =
        nova::fem::mark_dorfler_elements_3d(indicators, 0.8);
    NOVA_REQUIRE(eighty_percent.size() == 2U);
    NOVA_REQUIRE(eighty_percent[0] == 3U);
    NOVA_REQUIRE(eighty_percent[1] == 2U);

    const auto tied = nova::fem::mark_dorfler_elements_3d(
        std::vector<double>{2.0, 2.0, 1.0}, 0.5);
    NOVA_REQUIRE(tied.size() == 2U);
    NOVA_REQUIRE(tied[0] == 0U);
    NOVA_REQUIRE(tied[1] == 1U);

    NOVA_REQUIRE(nova::fem::mark_dorfler_elements_3d(
        std::vector<double>{0.0, 0.0, 0.0}, 0.5).empty());
    NOVA_REQUIRE(nova::fem::mark_dorfler_elements_3d(
        std::vector<double>{}, 0.5).empty());

    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::mark_dorfler_elements_3d(indicators, 0.0);
    }));
    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::mark_dorfler_elements_3d(indicators, 1.01);
    }));
    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::mark_dorfler_elements_3d(
            indicators, std::numeric_limits<double>::quiet_NaN());
    }));
    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::mark_dorfler_elements_3d(
            std::vector<double>{1.0, -1.0}, 0.5);
    }));
    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::mark_dorfler_elements_3d(
            std::vector<double>{1.0, std::numeric_limits<double>::infinity()}, 0.5);
    }));
}

NOVA_TEST(tetra_adaptive_refinement_is_local_conforming_and_transfers_fields) {
    const TetrahedralMesh3D coarse{
        {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {1.0, 1.0, 1.0}
        },
        {{0U, 1U, 2U, 3U}, {1U, 2U, 3U, 4U}}
    };
    coarse.validate();
    const auto coarse_boundary = nova::fem::extract_boundary_faces_3d(coarse);
    const auto refined = nova::fem::refine_tetrahedra_adaptively_3d(
        coarse, std::vector<std::size_t>{0U}, coarse_boundary);
    refined.mesh.validate();

    NOVA_REQUIRE(refined.coarse_node_count == coarse.nodes.size());
    NOVA_REQUIRE(refined.coarse_element_count == coarse.elements.size());
    NOVA_REQUIRE(refined.midpoint_parent_edges.size() == 6U);
    NOVA_REQUIRE(refined.mesh.nodes.size() == coarse.nodes.size() + 6U);
    NOVA_REQUIRE(refined.mesh.elements.size() == 12U);
    NOVA_REQUIRE(refined.parent_elements.size() == refined.mesh.elements.size());

    std::array<std::size_t, 2> child_counts{};
    for (const std::size_t parent : refined.parent_elements) {
        NOVA_REQUIRE(parent < child_counts.size());
        ++child_counts[parent];
    }
    NOVA_REQUIRE(child_counts[0] == 8U);
    NOVA_REQUIRE(child_counts[1] == 4U);

    const auto incidence = face_incidence(refined.mesh);
    for (const auto& [face, count] : incidence) {
        (void)face;
        NOVA_REQUIRE(count == 1U || count == 2U);
    }
    const auto midpoint = [&refined](std::size_t a, std::size_t b) {
        const std::array<std::size_t, 2> edge{{
            std::min(a, b), std::max(a, b)
        }};
        for (std::size_t i = 0; i < refined.midpoint_parent_edges.size(); ++i) {
            if (refined.midpoint_parent_edges[i] == edge) {
                return refined.coarse_node_count + i;
            }
        }
        throw std::runtime_error("Expected adaptive edge midpoint was not created.");
    };
    const std::size_t m12 = midpoint(1U, 2U);
    const std::size_t m23 = midpoint(2U, 3U);
    const std::size_t m31 = midpoint(1U, 3U);
    const std::set<std::size_t> shared_face_nodes{
        1U, 2U, 3U, m12, m23, m31
    };
    std::size_t shared_subface_count = 0U;
    for (const auto& [face, count] : incidence) {
        const bool lies_on_shared_face = std::all_of(
            face.begin(), face.end(),
            [&shared_face_nodes](const std::size_t node) {
                return shared_face_nodes.contains(node);
            });
        if (lies_on_shared_face) {
            ++shared_subface_count;
            NOVA_REQUIRE(count == 2U);
        }
    }
    NOVA_REQUIRE(shared_subface_count == 4U);

    const auto expected_boundary = nova::fem::extract_boundary_faces_3d(refined.mesh);
    NOVA_REQUIRE(expected_boundary.size() == refined.boundary_faces.size());
    std::map<FaceKey, BoundaryFace3D> expected_by_key;
    for (const auto& face : expected_boundary) {
        expected_by_key.emplace(canonical(face), face);
    }
    for (std::size_t i = 0; i < refined.boundary_faces.size(); ++i) {
        const auto found = expected_by_key.find(canonical(refined.boundary_faces[i]));
        NOVA_REQUIRE(found != expected_by_key.end());
        NOVA_REQUIRE(same_outward_winding(
            refined.boundary_faces[i], found->second));
        NOVA_REQUIRE(refined.parent_boundary_faces[i] < coarse_boundary.size());
    }

    const auto coarse_quality = nova::fem::evaluate_mesh_quality_3d(coarse);
    const auto refined_quality = nova::fem::evaluate_mesh_quality_3d(refined.mesh);
    NOVA_REQUIRE_NEAR(
        coarse_quality.total_volume, refined_quality.total_volume, 2.0e-13);
    NOVA_REQUIRE(refined_quality.minimum_mean_ratio > 0.0);

    std::vector<double> coarse_temperature;
    for (const auto& node : coarse.nodes) {
        coarse_temperature.push_back(
            2.0 + 3.0 * node.x - 4.0 * node.y + 5.0 * node.z);
    }
    const auto fine_temperature =
        nova::fem::prolongate_p1_refined_3d(refined, coarse_temperature);
    NOVA_REQUIRE(fine_temperature.size() == refined.mesh.nodes.size());
    for (std::size_t i = 0; i < fine_temperature.size(); ++i) {
        const auto& node = refined.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(
            fine_temperature[i],
            2.0 + 3.0 * node.x - 4.0 * node.y + 5.0 * node.z,
            2.0e-13);
    }

    const std::vector<double> coarse_cell_data{7.0, 11.0};
    const auto fine_cell_data =
        nova::fem::prolongate_cellwise_constant_refined_3d(
            refined, coarse_cell_data);
    NOVA_REQUIRE(fine_cell_data.size() == refined.parent_elements.size());
    long double coarse_integral = 0.0L;
    long double refined_integral = 0.0L;
    for (std::size_t i = 0; i < coarse_cell_data.size(); ++i) {
        coarse_integral += static_cast<long double>(coarse_cell_data[i]) *
            coarse_quality.elements[i].volume;
    }
    for (std::size_t i = 0; i < fine_cell_data.size(); ++i) {
        NOVA_REQUIRE(fine_cell_data[i] ==
            coarse_cell_data[refined.parent_elements[i]]);
        refined_integral += static_cast<long double>(fine_cell_data[i]) *
            refined_quality.elements[i].volume;
    }
    NOVA_REQUIRE_NEAR(
        static_cast<double>(coarse_integral),
        static_cast<double>(refined_integral), 2.0e-13);
}

NOVA_TEST(tetra_adaptive_refinement_handles_empty_and_invalid_marking) {
    const auto coarse = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    const auto unchanged = nova::fem::refine_tetrahedra_adaptively_3d(
        coarse.mesh, std::vector<std::size_t>{}, coarse.boundary_faces);
    NOVA_REQUIRE(unchanged.mesh.nodes.size() == coarse.mesh.nodes.size());
    NOVA_REQUIRE(unchanged.mesh.elements == coarse.mesh.elements);
    NOVA_REQUIRE(unchanged.parent_elements.size() == coarse.mesh.elements.size());
    NOVA_REQUIRE(unchanged.boundary_faces == coarse.boundary_faces);
    for (std::size_t i = 0; i < unchanged.parent_elements.size(); ++i) {
        NOVA_REQUIRE(unchanged.parent_elements[i] == i);
    }

    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::refine_tetrahedra_adaptively_3d(
            coarse.mesh, std::vector<std::size_t>{0U, 0U});
    }));
    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::refine_tetrahedra_adaptively_3d(
            coarse.mesh, std::vector<std::size_t>{coarse.mesh.elements.size()});
    }));
}

NOVA_TEST(tetra_adaptive_refinement_preserves_reversed_orientation) {
    TetrahedralMesh3D coarse{
        {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        },
        {{0U, 2U, 1U, 3U}}
    };
    coarse.validate();
    auto boundary = nova::fem::extract_boundary_faces_3d(coarse);
    for (auto& face : boundary) {
        std::rotate(face.begin(), face.begin() + 1, face.end());
    }
    const auto refined = nova::fem::refine_tetrahedra_adaptively_3d(
        coarse, std::vector<std::size_t>{0U}, boundary);
    refined.mesh.validate();
    NOVA_REQUIRE(refined.mesh.elements.size() == 8U);
    NOVA_REQUIRE(refined.boundary_faces.size() == 16U);
    NOVA_REQUIRE_NEAR(
        nova::fem::evaluate_mesh_quality_3d(refined.mesh).total_volume,
        1.0 / 6.0, 2.0e-14);
}


NOVA_TEST(tetra_adaptive_refinement_is_conforming_for_every_box_marking_subset) {
    const auto coarse = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    const auto coarse_quality =
        nova::fem::evaluate_mesh_quality_3d(coarse.mesh);
    const std::size_t subset_count =
        std::size_t{1U} << coarse.mesh.elements.size();

    for (std::size_t mask = 0; mask < subset_count; ++mask) {
        std::vector<std::size_t> marked;
        for (std::size_t element = 0; element < coarse.mesh.elements.size(); ++element) {
            if ((mask & (std::size_t{1U} << element)) != 0U) {
                marked.push_back(element);
            }
        }

        const auto refined = nova::fem::refine_tetrahedra_adaptively_3d(
            coarse.mesh, marked, coarse.boundary_faces);
        refined.mesh.validate();
        NOVA_REQUIRE(refined.parent_elements.size() == refined.mesh.elements.size());

        const auto incidence = face_incidence(refined.mesh);
        for (const auto& [face, count] : incidence) {
            (void)face;
            NOVA_REQUIRE(count == 1U || count == 2U);
        }

        const auto quality =
            nova::fem::evaluate_mesh_quality_3d(refined.mesh);
        NOVA_REQUIRE_NEAR(
            quality.total_volume, coarse_quality.total_volume, 3.0e-13);

        const auto expected_boundary =
            nova::fem::extract_boundary_faces_3d(refined.mesh);
        NOVA_REQUIRE(expected_boundary.size() == refined.boundary_faces.size());
        std::map<FaceKey, BoundaryFace3D> expected_by_key;
        for (const auto& face : expected_boundary) {
            expected_by_key.emplace(canonical(face), face);
        }
        for (std::size_t i = 0; i < refined.boundary_faces.size(); ++i) {
            const auto found = expected_by_key.find(canonical(refined.boundary_faces[i]));
            NOVA_REQUIRE(found != expected_by_key.end());
            NOVA_REQUIRE(same_outward_winding(
                refined.boundary_faces[i], found->second));
            NOVA_REQUIRE(refined.parent_boundary_faces[i] <
                coarse.boundary_faces.size());
        }
    }
}

NOVA_TEST(tetra_indicator_driven_refinement_matches_explicit_marking) {
    const auto coarse = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    NOVA_REQUIRE(coarse.mesh.elements.size() >= 2U);
    std::vector<double> indicators(coarse.mesh.elements.size(), 0.0);
    indicators[0] = 3.0;
    indicators[1] = 4.0;

    // The first ranked cell contributes 16/25 > theta=0.6.
    const auto result = nova::fem::refine_tetrahedra_by_indicators_3d(
        coarse.mesh, indicators, 0.6, coarse.boundary_faces);
    const std::vector<std::size_t> expected_marking{1U};
    NOVA_REQUIRE(result.marked_elements == expected_marking);
    const auto explicit_refinement = nova::fem::refine_tetrahedra_adaptively_3d(
        coarse.mesh, expected_marking, coarse.boundary_faces);
    NOVA_REQUIRE(
        result.refinement.mesh.elements == explicit_refinement.mesh.elements);
    NOVA_REQUIRE(
        result.refinement.parent_elements == explicit_refinement.parent_elements);
    NOVA_REQUIRE(
        result.refinement.midpoint_parent_edges ==
        explicit_refinement.midpoint_parent_edges);
    NOVA_REQUIRE(
        result.refinement.boundary_faces == explicit_refinement.boundary_faces);
    NOVA_REQUIRE(
        result.refinement.parent_boundary_faces ==
        explicit_refinement.parent_boundary_faces);

    const auto coarse_quality = nova::fem::evaluate_mesh_quality_3d(coarse.mesh);
    const auto fine_quality =
        nova::fem::evaluate_mesh_quality_3d(result.refinement.mesh);
    NOVA_REQUIRE_NEAR(
        fine_quality.total_volume, coarse_quality.total_volume, 3.0e-13);
}

NOVA_TEST(tetra_indicator_driven_refinement_zero_and_invalid_inputs) {
    const auto coarse = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    std::vector<double> indicators(coarse.mesh.elements.size(), 0.0);
    const auto unchanged = nova::fem::refine_tetrahedra_by_indicators_3d(
        coarse.mesh, indicators, 0.5, coarse.boundary_faces);
    NOVA_REQUIRE(unchanged.marked_elements.empty());
    NOVA_REQUIRE(unchanged.refinement.midpoint_parent_edges.empty());
    NOVA_REQUIRE(unchanged.refinement.mesh.elements == coarse.mesh.elements);
    NOVA_REQUIRE(unchanged.refinement.boundary_faces == coarse.boundary_faces);
    NOVA_REQUIRE(unchanged.refinement.parent_elements.size() ==
        coarse.mesh.elements.size());

    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::refine_tetrahedra_by_indicators_3d(
            coarse.mesh, std::vector<double>{1.0}, 0.5);
    }));
    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::refine_tetrahedra_by_indicators_3d(
            coarse.mesh, indicators, 0.0);
    }));
    indicators[0] = -1.0;
    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::refine_tetrahedra_by_indicators_3d(
            coarse.mesh, indicators, 0.5);
    }));
    indicators[0] = std::numeric_limits<double>::infinity();
    NOVA_REQUIRE(throws_mesh_error([&] {
        (void)nova::fem::refine_tetrahedra_by_indicators_3d(
            coarse.mesh, indicators, 0.5);
    }));
}
