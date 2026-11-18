#include "test_framework.hpp"

#include "nova/io/vtk_legacy_2d.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

nova::fem::TriangleMesh2D two_triangle_mesh() {
    nova::fem::TriangleMesh2D mesh;
    mesh.nodes = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    mesh.elements = {{0U, 1U, 2U}, {0U, 2U, 3U}};
    return mesh;
}

void require_scalar_values(
    const std::string& vtk,
    const std::vector<double>& expected
) {
    const std::string header = "LOOKUP_TABLE default\n";
    const auto position = vtk.find(header);
    NOVA_REQUIRE(position != std::string::npos);
    std::istringstream data(vtk.substr(position + header.size()));
    for (const double value : expected) {
        double actual = 0.0;
        NOVA_REQUIRE(static_cast<bool>(data >> actual));
        NOVA_REQUIRE_NEAR(actual, value, 1.0e-13);
    }
    double extra = 0.0;
    NOVA_REQUIRE(!(data >> extra));
}

} // namespace

NOVA_TEST(vtk_export_p1_triangle_connectivity_and_scalar_order) {
    const auto mesh = two_triangle_mesh();
    const nova::fem::H1FunctionSpace2D space(
        mesh, nova::fem::H1TriangleOrder::P1
    );
    const auto values = space.interpolate(
        [](double x, double y) { return 2.0 * x + y; }
    );
    std::ostringstream out;
    nova::io::write_legacy_vtk_scalar_2d(out, space, values);
    const auto vtk = out.str();
    NOVA_REQUIRE(vtk.find("# vtk DataFile Version 3.0\n") == 0U);
    NOVA_REQUIRE(vtk.find("DATASET UNSTRUCTURED_GRID\n") != std::string::npos);
    NOVA_REQUIRE(vtk.find("POINTS 4 double\n") != std::string::npos);
    NOVA_REQUIRE(vtk.find(
        "CELLS 2 8\n3 0 1 2\n3 0 2 3\n"
    ) != std::string::npos);
    NOVA_REQUIRE(vtk.find(
        "CELL_TYPES 2\n5\n5\n"
    ) != std::string::npos);
    NOVA_REQUIRE(vtk.find(
        "POINT_DATA 4\nSCALARS temperature double 1\n"
    ) != std::string::npos);
    require_scalar_values(vtk, values);
}

NOVA_TEST(vtk_export_p2_quadratic_cells_share_midpoint_dofs) {
    const auto mesh = two_triangle_mesh();
    const nova::fem::H1FunctionSpace2D space(
        mesh, nova::fem::H1TriangleOrder::P2
    );
    const auto values = space.interpolate(
        [](double x, double y) { return x * x + y; }
    );
    NOVA_REQUIRE(values.size() == 9U);
    std::ostringstream out;
    nova::io::write_legacy_vtk_scalar_2d(
        out, space, values, "temperature_K"
    );
    const auto vtk = out.str();
    NOVA_REQUIRE(vtk.find("POINTS 9 double\n") != std::string::npos);
    // Global P2 point IDs: 0-3 vertices, 4=(0,1), 5=(1,2),
    // 6=(2,0), 7=(2,3), 8=(3,0). ID 6 is shared by both cells.
    NOVA_REQUIRE(vtk.find(
        "CELLS 2 14\n6 0 1 2 4 5 6\n6 0 2 3 6 7 8\n"
    ) != std::string::npos);
    NOVA_REQUIRE(vtk.find(
        "CELL_TYPES 2\n22\n22\n"
    ) != std::string::npos);
    NOVA_REQUIRE(vtk.find(
        "SCALARS temperature_K double 1\n"
    ) != std::string::npos);
    require_scalar_values(vtk, values);
}

NOVA_TEST(vtk_export_rejects_mismatched_and_nonfinite_scalar_values) {
    const auto mesh = two_triangle_mesh();
    const nova::fem::H1FunctionSpace2D space(
        mesh, nova::fem::H1TriangleOrder::P2
    );
    auto values = space.interpolate(
        [](double, double) { return 1.0; }
    );
    std::ostringstream out;
    values.pop_back();
    bool wrong_count = false;
    try {
        nova::io::write_legacy_vtk_scalar_2d(out, space, values);
    } catch (const std::invalid_argument&) {
        wrong_count = true;
    }
    NOVA_REQUIRE(wrong_count);
    NOVA_REQUIRE(out.str().empty());

    values.push_back(std::numeric_limits<double>::quiet_NaN());
    bool invalid_value = false;
    try {
        nova::io::write_legacy_vtk_scalar_2d(out, space, values);
    } catch (const std::invalid_argument&) {
        invalid_value = true;
    }
    NOVA_REQUIRE(invalid_value);
    NOVA_REQUIRE(out.str().empty());
}

NOVA_TEST(vtk_export_rejects_invalid_field_name_and_unwritable_stream) {
    const auto mesh = two_triangle_mesh();
    const nova::fem::H1FunctionSpace2D space(
        mesh, nova::fem::H1TriangleOrder::P1
    );
    const auto values = space.interpolate(
        [](double, double) { return 1.0; }
    );
    std::ostringstream out;
    bool invalid_name = false;
    try {
        nova::io::write_legacy_vtk_scalar_2d(
            out, space, values, "temperature\nINJECTED"
        );
    } catch (const std::invalid_argument&) {
        invalid_name = true;
    }
    NOVA_REQUIRE(invalid_name);
    NOVA_REQUIRE(out.str().empty());

    out.setstate(std::ios::badbit);
    bool bad_stream = false;
    try {
        nova::io::write_legacy_vtk_scalar_2d(out, space, values);
    } catch (const std::runtime_error&) {
        bad_stream = true;
    }
    NOVA_REQUIRE(bad_stream);
}
