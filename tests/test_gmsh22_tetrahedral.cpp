#include "test_framework.hpp"

#include "nova/io/gmsh22_tetrahedral.hpp"
#include "nova/fem/diffusion_3d.hpp"

#include <algorithm>
#include <array>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

[[nodiscard]] std::string sample_mesh() {
    return R"msh($MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
2 7 "wall"
3 42 "solid"
$EndPhysicalNames
$Nodes
4
90 0 0 1
10 0 0 0
30 1 0 0
50 0 1 0
$EndNodes
$Elements
4
101 2 2 7 1 10 30 50
200 4 2 42 1 10 30 50 90
103 1 2 19 1 10 30
104 15 2 99 1 10
$EndElements
)msh";
}

[[nodiscard]] nova::io::GmshTetrahedralMesh3D parse(const std::string& content) {
    std::istringstream stream(content);
    return nova::io::read_gmsh22_tetrahedra(stream);
}

void reject(const std::string& content) {
    bool rejected = false;
    try {
        static_cast<void>(parse(content));
    } catch (const std::invalid_argument&) {
        rejected = true;
    } catch (const std::out_of_range&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}

[[nodiscard]] std::string replaced(std::string source,
                                   const std::string& needle,
                                   const std::string& replacement) {
    const auto at = source.find(needle);
    NOVA_REQUIRE(at != std::string::npos);
    source.replace(at, needle.size(), replacement);
    return source;
}

} // namespace

NOVA_TEST(gmsh22_import_remaps_sparse_node_ids_and_preserves_physical_tags) {
    const auto imported = parse(sample_mesh());
    imported.mesh.validate();
    NOVA_REQUIRE(imported.mesh.nodes.size() == 4U);
    NOVA_REQUIRE(imported.mesh.elements.size() == 1U);
    NOVA_REQUIRE(imported.mesh.elements[0] == (nova::fem::Tetrahedron3D{1, 2, 3, 0}));
    NOVA_REQUIRE(imported.volume_physical_tags == std::vector<int>{42});
    NOVA_REQUIRE(imported.exterior_faces.size() == 4U);
    NOVA_REQUIRE(imported.exterior_physical_tags.size() == 4U);
    NOVA_REQUIRE(std::count(imported.exterior_physical_tags.begin(),
                            imported.exterior_physical_tags.end(), 7) == 1);
    NOVA_REQUIRE(std::count(imported.exterior_physical_tags.begin(),
                            imported.exterior_physical_tags.end(), 0) == 3);
    for (std::size_t i = 0; i < imported.exterior_faces.size(); ++i) {
        const auto& face = imported.exterior_faces[i];
        const bool on_bottom = std::all_of(face.begin(), face.end(), [&](std::size_t node) {
            return imported.mesh.nodes[node].z == 0.0;
        });
        NOVA_REQUIRE((imported.exterior_physical_tags[i] == 7) == on_bottom);
    }
}

NOVA_TEST(gmsh22_import_untagged_boundary_faces_support_diffusion_solver) {
    std::string content = replaced(
        sample_mesh(), "$Elements\n4\n", "$Elements\n3\n"
    );
    content = replaced(content, "101 2 2 7 1 10 30 50\n", "");
    const auto imported = parse(content);
    NOVA_REQUIRE(std::all_of(
        imported.exterior_physical_tags.begin(),
        imported.exterior_physical_tags.end(), [](int tag) { return tag == 0; }
    ));

    nova::fem::ScalarDiffusion3DProblem problem;
    problem.mesh = imported.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    for (const auto& face : imported.exterior_faces) {
        problem.dirichlet_faces.push_back({
            face,
            [](double x, double y, double z) { return x + 2.0 * y - 3.0 * z; }
        });
    }
    const auto result = nova::fem::solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
        const auto& p = problem.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(result.values[i], p.x + 2.0 * p.y - 3.0 * p.z, 1.0e-12);
    }
}

NOVA_TEST(gmsh22_import_rejects_unsupported_versions_types_and_bad_sections) {
    reject(replaced(sample_mesh(), "2.2 0 8", "4.1 0 8"));
    reject(replaced(sample_mesh(), "2.2 0 8", "2.2 1 8"));
    reject(replaced(sample_mesh(), "2.2 0 8", "2.2 0 4"));
    reject(replaced(sample_mesh(), "200 4 2", "200 11 2"));
    reject(replaced(sample_mesh(), "$EndNodes", "$EndNodex"));
    reject(replaced(sample_mesh(), "$EndElements", "$EndElementx"));
    reject(sample_mesh() + "surprise\n");
}

NOVA_TEST(gmsh22_import_rejects_corrupt_nodes_elements_and_faces) {
    reject(replaced(sample_mesh(), "50 0 1 0", "30 0 1 0")); // duplicate node ID
    reject(replaced(sample_mesh(), "200 4 2 42 1 10", "200 4 2 42 1 999"));
    reject(replaced(sample_mesh(), "200 4 2 42 1", "200 4 2 42 1 500"));
    reject(replaced(sample_mesh(), "101 2 2 7 1 10 30 50",
                   "200 2 2 7 1 10 30 50")); // repeated element ID
    reject(replaced(sample_mesh(), "101 2 2 7 1 10 30 50",
                   "101 2 2 7 1 10 10 50")); // invalid boundary triangle
    reject(replaced(sample_mesh(), "90 0 0 1", "90 0 0 nan"));
    reject(replaced(sample_mesh(), "90 0 0 1", "90 0 0 0")); // zero volume
    const auto repeated = replaced(
        replaced(sample_mesh(), "$Elements\n4\n", "$Elements\n5\n"),
        "$EndElements", "105 2 2 7 1 50 30 10\n$EndElements"
    );
    reject(repeated);
}

NOVA_TEST(gmsh22_import_rejects_interior_surface_triangles) {
    const std::string file = R"msh($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
5
10 0 0 0
20 1 0 0
30 0 1 0
40 0 0 1
50 0 0 -1
$EndNodes
$Elements
3
1 4 0 10 20 30 40
2 4 0 10 30 20 50
3 2 1 7 10 20 30
$EndElements
)msh";
    reject(file);
}

NOVA_TEST(gmsh22_import_rejects_unused_3d_vertices) {
    std::string file = replaced(
        sample_mesh(), "$Nodes\n4\n", "$Nodes\n5\n"
    );
    file = replaced(file, "$EndNodes", "100 0.5 0.5 0.5\n$EndNodes");
    reject(file);
}
