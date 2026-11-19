#include "test_framework.hpp"

#include "nova/io/vtk_legacy_3d.hpp"

#include <cmath>
#include <ios>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

nova::fem::TetrahedralMesh3D one_tetrahedron() {
    return {
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
         {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}},
        {{0U, 1U, 2U, 3U}}
    };
}

void require_scalar_values(const std::string& vtk,
                           const std::vector<double>& expected) {
    const std::string marker = "LOOKUP_TABLE default\n";
    const auto location = vtk.find(marker);
    NOVA_REQUIRE(location != std::string::npos);
    std::istringstream values(vtk.substr(location + marker.size()));
    for (const double expected_value : expected) {
        double actual = 0.0;
        NOVA_REQUIRE(static_cast<bool>(values >> actual));
        NOVA_REQUIRE_NEAR(actual, expected_value, 1.0e-14);
    }
    double extra = 0.0;
    NOVA_REQUIRE(!(values >> extra));
}

} // namespace

NOVA_TEST(vtk_3d_p1_exports_tetrahedral_connectivity_and_scalar_order) {
    const auto mesh = one_tetrahedron();
    const nova::fem::H1FunctionSpace3D space(mesh);
    const auto values = space.interpolate(
        [](double x, double y, double z) { return 2*x - y + 3*z; }
    );
    std::ostringstream output;
    nova::io::write_legacy_vtk_scalar_3d(output, space, values);
    const auto vtk = output.str();
    NOVA_REQUIRE(vtk.find("# vtk DataFile Version 3.0\n") == 0U);
    NOVA_REQUIRE(vtk.find("DATASET UNSTRUCTURED_GRID\n") != std::string::npos);
    NOVA_REQUIRE(vtk.find("POINTS 4 double\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n")
                 != std::string::npos);
    NOVA_REQUIRE(vtk.find("CELLS 1 5\n4 0 1 2 3\n") != std::string::npos);
    NOVA_REQUIRE(vtk.find("CELL_TYPES 1\n10\n") != std::string::npos);
    NOVA_REQUIRE(vtk.find("POINT_DATA 4\nSCALARS temperature double 1\n")
                 != std::string::npos);
    require_scalar_values(vtk, values);
}

NOVA_TEST(vtk_3d_p1_preserves_global_vertex_ids_over_shared_tetrahedra) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U
    );
    const nova::fem::H1FunctionSpace3D space(box.mesh);
    const auto values = space.interpolate(
        [](double x, double y, double z) {
            return 0.125 + x + 0.5*y - 0.25*z;
        }
    );
    std::ostringstream output;
    nova::io::write_legacy_vtk_scalar_3d(output, space, values, "temperature_K");
    const auto vtk = output.str();
    NOVA_REQUIRE(vtk.find("POINTS 8 double\n") != std::string::npos);
    NOVA_REQUIRE(vtk.find("CELL_TYPES 6\n10\n10\n10\n10\n10\n10\n")
                 != std::string::npos);
    NOVA_REQUIRE(vtk.find("SCALARS temperature_K double 1\n")
                 != std::string::npos);

    const std::string marker = "CELLS 6 30\n";
    const auto location = vtk.find(marker);
    NOVA_REQUIRE(location != std::string::npos);
    std::istringstream cells(vtk.substr(location + marker.size()));
    for (std::size_t e = 0; e < box.mesh.elements.size(); ++e) {
        std::size_t count = 0;
        NOVA_REQUIRE(static_cast<bool>(cells >> count));
        NOVA_REQUIRE(count == 4U);
        for (std::size_t j = 0; j < 4U; ++j) {
            std::size_t node = 0;
            NOVA_REQUIRE(static_cast<bool>(cells >> node));
            NOVA_REQUIRE(node == space.dofs().element_dofs(e)[j]);
        }
    }
    require_scalar_values(vtk, values);
}

NOVA_TEST(vtk_3d_p1_rejects_bad_field_data_without_writing) {
    const auto mesh = one_tetrahedron();
    const nova::fem::H1FunctionSpace3D space(mesh);
    std::ostringstream output;
    std::vector<double> values{0.0, 1.0, 2.0};
    bool rejected = false;
    try {
        nova::io::write_legacy_vtk_scalar_3d(output, space, values);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
    NOVA_REQUIRE(output.str().empty());

    values.push_back(std::numeric_limits<double>::infinity());
    rejected = false;
    try {
        nova::io::write_legacy_vtk_scalar_3d(output, space, values);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
    NOVA_REQUIRE(output.str().empty());
}

NOVA_TEST(vtk_3d_p1_rejects_invalid_name_and_failed_stream) {
    const auto mesh = one_tetrahedron();
    const nova::fem::H1FunctionSpace3D space(mesh);
    const auto values = space.interpolate(
        [](double, double, double) { return 1.0; }
    );
    std::ostringstream output;
    bool rejected = false;
    try {
        nova::io::write_legacy_vtk_scalar_3d(
            output, space, values, "T\nCELL_DATA 10"
        );
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
    NOVA_REQUIRE(output.str().empty());

    output.setstate(std::ios::badbit);
    bool failed_stream = false;
    try {
        nova::io::write_legacy_vtk_scalar_3d(output, space, values);
    } catch (const std::runtime_error&) {
        failed_stream = true;
    }
    NOVA_REQUIRE(failed_stream);
}
