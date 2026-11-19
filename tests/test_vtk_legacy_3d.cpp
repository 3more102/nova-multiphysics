#include "test_framework.hpp"
#include "nova/io/vtk_legacy_3d.hpp"

#include <array>
#include <cmath>
#include <ios>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

nova::fem::TetrahedralMesh3D one_tetrahedron() {
    return {{
        {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}
    }, {{0,1,2,3}}};
}

template <typename Exception, typename F>
bool throws_as(F&& action) {
    try { action(); } catch (const Exception&) { return true; }
    return false;
}

} // namespace

NOVA_TEST(vtk3d_p1_tetrahedral_connectivity_scalar_order_and_precision) {
    const auto mesh = one_tetrahedron();
    const std::vector values{1.2345678901234567, -0.25, 0.0, 2.5};
    std::ostringstream output;
    output.precision(3);
    nova::io::write_legacy_vtk_scalar_3d(output, mesh, values, "temperature_K");
    NOVA_REQUIRE(output.precision() == 3);
    const std::string result = output.str();
    NOVA_REQUIRE(result.find("# vtk DataFile Version 3.0\n") == 0U);
    NOVA_REQUIRE(result.find("DATASET UNSTRUCTURED_GRID\n") != std::string::npos);
    NOVA_REQUIRE(result.find("POINTS 4 double\n") != std::string::npos);
    NOVA_REQUIRE(result.find("CELLS 1 5\n4 0 1 2 3\n") != std::string::npos);
    NOVA_REQUIRE(result.find("CELL_TYPES 1\n10\n") != std::string::npos);
    NOVA_REQUIRE(result.find(
        "POINT_DATA 4\nSCALARS temperature_K double 1\nLOOKUP_TABLE default\n"
    ) != std::string::npos);
    const std::string marker = "LOOKUP_TABLE default\n";
    const auto position = result.find(marker);
    NOVA_REQUIRE(position != std::string::npos);
    std::istringstream read(result.substr(position + marker.size()));
    for (double expected : values) {
        double parsed = 0.0;
        NOVA_REQUIRE(static_cast<bool>(read >> parsed));
        NOVA_REQUIRE(parsed == expected);
    }
    double extra = 0.0;
    NOVA_REQUIRE(!(read >> extra));
}

NOVA_TEST(vtk3d_two_cells_keep_shared_global_vertex_ids) {
    // Cells share face (1,2,3), with their opposing vertices on opposite sides.
    nova::fem::TetrahedralMesh3D mesh = {{
        {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,1,1}
    }, {{0,1,2,3}, {1,2,3,4}}};
    std::ostringstream output;
    nova::io::write_legacy_vtk_scalar_3d(
        output, mesh, std::vector{0., 1., 2., 3., 4.});
    NOVA_REQUIRE(output.str().find(
        "CELLS 2 10\n4 0 1 2 3\n4 1 2 3 4\n") != std::string::npos);
    NOVA_REQUIRE(output.str().find(
        "CELL_TYPES 2\n10\n10\n") != std::string::npos);
}

NOVA_TEST(vtk3d_rejects_nonfinite_mismatched_or_injected_scalar_data) {
    const auto mesh = one_tetrahedron();
    std::ostringstream output;
    NOVA_REQUIRE(throws_as<std::invalid_argument>([&] {
        nova::io::write_legacy_vtk_scalar_3d(output, mesh, {1.,2.,3.});
    }));
    NOVA_REQUIRE(output.str().empty());
    NOVA_REQUIRE(throws_as<std::invalid_argument>([&] {
        nova::io::write_legacy_vtk_scalar_3d(output, mesh,
            {1.,2.,3.,std::numeric_limits<double>::quiet_NaN()});
    }));
    NOVA_REQUIRE(output.str().empty());
    NOVA_REQUIRE(throws_as<std::invalid_argument>([&] {
        nova::io::write_legacy_vtk_scalar_3d(output, mesh,
            {1.,2.,3.,4.}, "heat\nCELL_TYPES 0");
    }));
    NOVA_REQUIRE(output.str().empty());
}

NOVA_TEST(vtk3d_rejects_invalid_mesh_and_unwritable_stream) {
    auto mesh = one_tetrahedron();
    mesh.elements.front()[1] = 0U;
    std::ostringstream output;
    NOVA_REQUIRE(throws_as<std::exception>([&] {
        nova::io::write_legacy_vtk_scalar_3d(output, mesh, {1.,2.,3.,4.});
    }));
    NOVA_REQUIRE(output.str().empty());

    output.setstate(std::ios::badbit);
    NOVA_REQUIRE(throws_as<std::runtime_error>([&] {
        nova::io::write_legacy_vtk_scalar_3d(output, one_tetrahedron(), {1.,2.,3.,4.});
    }));
}

NOVA_TEST(vtk3d_tagged_conduction_validates_cell_data_before_emitting_geometry) {
    const auto mesh = one_tetrahedron();
    const std::vector<double> temperature(4, 300.0);
    const std::array tags{12};
    std::array<nova::fem::P1TetrahedralFlux3D, 1> fields{};
    std::ostringstream output;
    NOVA_REQUIRE(throws_as<std::invalid_argument>([&] {
        nova::io::write_legacy_vtk_conduction_3d(output, mesh, temperature, {}, fields);
    }));
    NOVA_REQUIRE(throws_as<std::invalid_argument>([&] {
        nova::io::write_legacy_vtk_conduction_3d(output, mesh, temperature, tags, {});
    }));
    fields[0].flux.z = std::numeric_limits<double>::infinity();
    NOVA_REQUIRE(throws_as<std::invalid_argument>([&] {
        nova::io::write_legacy_vtk_conduction_3d(output, mesh, temperature, tags, fields);
    }));
    NOVA_REQUIRE(output.str().empty());

    fields[0].flux.z = 1e-50;
    output << std::hex << std::fixed;
    const auto flags = output.flags();
    nova::io::write_legacy_vtk_conduction_3d(output, mesh, temperature, tags, fields);
    NOVA_REQUIRE(output.flags() == flags);
    NOVA_REQUIRE(output.str().find("LOOKUP_TABLE default\n12\nVECTORS heat_flux double\n") != std::string::npos);
    NOVA_REQUIRE(output.str().find("e-50") != std::string::npos);
}


NOVA_TEST(vtk3d_heat_generation_field_is_precise_optional_and_prevalidated) {
    const auto mesh = one_tetrahedron();
    const std::vector<double> temperature(4, 300.0);
    const std::array tags{11};
    const std::array fields{nova::fem::P1TetrahedralFlux3D{}};
    const std::array sources{120.125};
    std::ostringstream output;
    const std::array<double, 2> wrong_size{120.0, 0.0};
    NOVA_REQUIRE(throws_as<std::invalid_argument>([&] {
        nova::io::write_legacy_vtk_conduction_3d(
            output, mesh, temperature, tags, fields, wrong_size);
    }));
    NOVA_REQUIRE(output.str().empty());
    const std::array nonfinite{std::numeric_limits<double>::quiet_NaN()};
    NOVA_REQUIRE(throws_as<std::invalid_argument>([&] {
        nova::io::write_legacy_vtk_conduction_3d(
            output, mesh, temperature, tags, fields, nonfinite);
    }));
    NOVA_REQUIRE(output.str().empty());
    output.precision(4);
    nova::io::write_legacy_vtk_conduction_3d(
        output, mesh, temperature, tags, fields, sources);
    NOVA_REQUIRE(output.precision() == 4);
    NOVA_REQUIRE(output.str().find(
        "SCALARS material_id int 1\nLOOKUP_TABLE default\n11\n"
        "SCALARS heat_generation double 1\nLOOKUP_TABLE default\n"
        "120.125\nVECTORS heat_flux double\n") != std::string::npos);
    std::ostringstream legacy;
    nova::io::write_legacy_vtk_conduction_3d(
        legacy, mesh, temperature, tags, fields);
    NOVA_REQUIRE(legacy.str().find("heat_generation") == std::string::npos);
}
