#include "test_framework.hpp"

#include "nova/io/gmsh_tetrahedral.hpp"
#include "nova/fem/diffusion_3d.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
const std::string fixture = R"msh($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
2 101 "fixed-temperature"
3 11 "solid"
$EndPhysicalNames
$Entities
0 0 1 1
7 0 0 0 1 1 0 1 101 0
9 0 0 0 1 1 1 1 11 0
$EndEntities
$Nodes
2 4 10 90
2 7 0 3
10
30
50
0 0 0
1 0 0
0 1 0
3 9 0 1
90
0 0 1
$EndNodes
$Elements
2 2 101 200
2 7 2 1
101 10 30 50
3 9 4 1
200 10 30 50 90
$EndElements
)msh";

nova::io::GmshTetrahedralMesh3D parse(const std::string& input) {
    std::istringstream stream(input);
    return nova::io::read_gmsh_tetrahedra(stream);
}

std::string replace(std::string input, const std::string& from, const std::string& to) {
    const auto at = input.find(from);
    NOVA_REQUIRE(at != std::string::npos);
    input.replace(at, from.size(), to);
    return input;
}

void reject(const std::string& input) {
    bool rejected = false;
    try { static_cast<void>(parse(input)); }
    catch (const std::invalid_argument&) { rejected = true; }
    catch (const std::out_of_range&) { rejected = true; }
    NOVA_REQUIRE(rejected);
}
} // namespace

NOVA_TEST(gmsh41_preserves_sparse_nodes_entity_physical_tags_and_outward_faces) {
    const auto mesh = parse(fixture);
    mesh.mesh.validate();
    NOVA_REQUIRE(mesh.mesh.nodes.size() == 4);
    NOVA_REQUIRE(mesh.mesh.elements.size() == 1);
    NOVA_REQUIRE(mesh.mesh.elements[0] == (nova::fem::Tetrahedron3D{0,1,2,3}));
    NOVA_REQUIRE(mesh.volume_physical_tags.size() == 1);
    NOVA_REQUIRE(mesh.volume_physical_tags[0] == 11);
    NOVA_REQUIRE(mesh.physical_names.at({2,101}) == "fixed-temperature");
    NOVA_REQUIRE(mesh.physical_names.at({3,11}) == "solid");
    NOVA_REQUIRE(mesh.exterior_physical_tags.size() == 4);
    NOVA_REQUIRE(std::count(mesh.exterior_physical_tags.begin(),
                            mesh.exterior_physical_tags.end(), 101) == 1);
    for (std::size_t i = 0; i < mesh.exterior_faces.size(); ++i) {
        const auto face = mesh.exterior_faces[i];
        const bool bottom = std::all_of(face.begin(), face.end(), [&](std::size_t n) {
            return mesh.mesh.nodes[n].z == 0.0;
        });
        NOVA_REQUIRE((mesh.exterior_physical_tags[i] == 101) == bottom);
    }
}

NOVA_TEST(gmsh41_thermal_solve_with_named_boundary_and_material) {
    const auto imported = parse(fixture);
    nova::fem::ScalarDiffusion3DProblem p;
    p.mesh = imported.mesh;
    p.cell_diffusivity = {2.0};
    p.source = [](double, double, double) { return 0.0; };
    for (std::size_t i = 0; i < imported.exterior_faces.size(); ++i) {
        if (imported.exterior_physical_tags[i] == 101) {
            p.dirichlet_faces.push_back({
                imported.exterior_faces[i],
                [](double, double, double) { return 320.0; }
            });
        }
    }
    const auto result = nova::fem::solve_scalar_diffusion_3d(p);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (double temperature : result.values) {
        NOVA_REQUIRE_NEAR(temperature, 320.0, 1e-9);
    }
}

NOVA_TEST(gmsh41_missing_entities_and_parametric_nodes_are_supported) {
    std::string plain = replace(fixture,
        "$Entities\n0 0 1 1\n7 0 0 0 1 1 0 1 101 0\n9 0 0 0 1 1 1 1 11 0\n$EndEntities\n", "");
    const auto untagged = parse(plain);
    NOVA_REQUIRE(untagged.volume_physical_tags[0] == 0);
    NOVA_REQUIRE(std::all_of(untagged.exterior_physical_tags.begin(),
                             untagged.exterior_physical_tags.end(),
                             [](int tag) { return tag == 0; }));
    auto parametric = replace(fixture, "2 7 0 3\n", "2 7 1 3\n");
    parametric = replace(parametric, "0 0 0\n1 0 0\n0 1 0\n",
                        "0 0 0 0 0\n1 0 0 1 0\n0 1 0 0 1\n");
    const auto restored = parse(parametric);
    NOVA_REQUIRE(restored.mesh.nodes.size() == 4);
}

NOVA_TEST(gmsh41_rejects_invalid_counts_missing_nodes_and_unsupported_elements) {
    reject(replace(fixture, "2 4 10 90", "2 5 10 90"));
    reject(replace(fixture, "2 4 10 90", "2 4 10 91"));
    reject(replace(fixture, "200 10 30 50 90", "200 10 30 50 999"));
    reject(replace(fixture, "3 9 4 1", "3 9 11 1"));
    reject(replace(fixture, "4.1 0 8", "4.1 1 8"));
    reject(replace(fixture, "9 0 0 0 1 1 1 1 11 0",
                            "9 0 0 0 1 1 1 2 11 12 0"));
    reject(replace(fixture, "2 7 2 1", "3 7 2 1"));
    reject(replace(fixture, "$EndElements\n", "$EndElements\n$Periodic\n"));
    reject(replace(fixture, "3 9 0 1", "3 100 0 1"));
}

NOVA_TEST(gmsh41_rejects_duplicate_sparse_identifiers) {
    reject(replace(fixture, "30\n50", "10\n50"));
    reject(replace(fixture, "200 10 30 50 90", "101 10 30 50 90"));
}

NOVA_TEST(gmsh41_generic_import_does_not_regress_msh22) {
    const std::string msh22 = R"msh($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
90 0 0 1
10 0 0 0
30 1 0 0
50 0 1 0
$EndNodes
$Elements
1
200 4 2 11 9 10 30 50 90
$EndElements
)msh";
    const auto mesh = parse(msh22);
    NOVA_REQUIRE(mesh.mesh.nodes.size() == 4);
    NOVA_REQUIRE(mesh.mesh.elements.size() == 1);
    NOVA_REQUIRE(mesh.volume_physical_tags[0] == 11);
}
