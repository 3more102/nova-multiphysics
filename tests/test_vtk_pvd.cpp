#include "test_framework.hpp"
#include "nova/io/vtk_pvd.hpp"

#include <array>
#include <ios>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

NOVA_TEST(vtk_pvd_collection_writes_valid_timed_datasets_and_escapes_names) {
    const std::array<nova::io::VtkTimeSeriesEntry, 3> entries{{
        {-2.0, "frames/initial.vtk"},
        {0.125, "frames/step_0001.vtk"},
        {0.25, "frames/end & finale's.vtu"}
    }};
    std::ostringstream output;
    output.precision(3);
    nova::io::write_vtk_pvd_collection(output, entries);
    NOVA_REQUIRE(output.precision() == 3);
    const auto xml = output.str();
    NOVA_REQUIRE(xml.find("<?xml version=\"1.0\"?>\n") == 0U);
    NOVA_REQUIRE(xml.find("<VTKFile type=\"Collection\"") != std::string::npos);
    NOVA_REQUIRE(xml.find("timestep=\"-2\" group=\"\" part=\"0\" file=\"frames/initial.vtk\"") != std::string::npos);
    NOVA_REQUIRE(xml.find("timestep=\"0.125\" group=\"\" part=\"0\" file=\"frames/step_0001.vtk\"") != std::string::npos);
    NOVA_REQUIRE(xml.find("timestep=\"0.25\" group=\"\" part=\"0\" file=\"frames/end &amp; finale&apos;s.vtu\"") != std::string::npos);
    NOVA_REQUIRE(xml.find("  </Collection>\n</VTKFile>\n") != std::string::npos);
}

NOVA_TEST(vtk_pvd_collection_prevalidates_strict_time_and_safe_file_paths) {
    const auto rejects = [](std::vector<nova::io::VtkTimeSeriesEntry> entries) {
        std::ostringstream output;
        try {
            nova::io::write_vtk_pvd_collection(output, entries);
        } catch (const std::invalid_argument&) {
            return output.str().empty();
        }
        return false;
    };

    NOVA_REQUIRE(rejects({}));
    NOVA_REQUIRE(rejects({{0.0, "one.vtk"}, {0.0, "two.vtk"}}));
    NOVA_REQUIRE(rejects({{0.0, "one.vtk"}, {-1.0, "two.vtk"}}));
    NOVA_REQUIRE(rejects({{0.0, "one.vtk"}, {1.0, "one.vtk"}}));
    NOVA_REQUIRE(rejects({{std::numeric_limits<double>::infinity(), "one.vtk"}}));
    NOVA_REQUIRE(rejects({{std::numeric_limits<double>::quiet_NaN(), "one.vtk"}}));
    for (const auto* bad : {
             "../step.vtk", "./step.vtk", "/root/step.vtk", "frames//step.vtk",
             "frames/../step.vtk", "C:/step.vtk", "frames\\step.vtk",
             "frames/step.vtk\n", "frames/step.txt", "frames/step?.vtk",
             "frames/step<evil>.vtk", "frames/"
         }) {
        NOVA_REQUIRE(rejects({{0.0, bad}}));
    }
}

NOVA_TEST(vtk_pvd_collection_rejects_bad_output_stream_before_writing) {
    std::ostringstream output;
    output.setstate(std::ios::badbit);
    const std::array<nova::io::VtkTimeSeriesEntry, 1> entries{{
        {0.0, "frame_0000.vtk"}
    }};
    bool caught = false;
    try {
        nova::io::write_vtk_pvd_collection(output, entries);
    } catch (const std::runtime_error&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}
