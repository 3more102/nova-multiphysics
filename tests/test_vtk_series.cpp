#include "test_framework.hpp"
#include "nova/io/vtk_series.hpp"

#include <cmath>
#include <iomanip>
#include <ios>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
template <typename Exception, typename Fn>
bool raises(Fn&& function) {
    try { function(); } catch (const Exception&) { return true; }
    return false;
}
} // namespace

NOVA_TEST(pvd_collection_preserves_time_precision_and_stream_format) {
    std::ostringstream output;
    output << std::fixed << std::setprecision(3);
    const auto original_flags = output.flags();
    nova::io::write_pvd_collection(output, {
        {0.0, "cooling_000000.vtk"},
        {0.12345678901234566, "cooling_000001.vtk"},
        {0.25, "cooling_000002.vtk"}
    });
    NOVA_REQUIRE(output.flags() == original_flags);
    NOVA_REQUIRE(output.precision() == 3);
    const auto text = output.str();
    NOVA_REQUIRE(text.find("<?xml version=\"1.0\"?>\n") == 0U);
    NOVA_REQUIRE(text.find("type=\"Collection\"") != std::string::npos);
    NOVA_REQUIRE(text.find("timestep=\"0\"") != std::string::npos);
    NOVA_REQUIRE(text.find("timestep=\"0.12345678901234566\"") != std::string::npos);
    NOVA_REQUIRE(text.find("file=\"cooling_000002.vtk\"") != std::string::npos);
    NOVA_REQUIRE(text.find("</Collection>\n</VTKFile>\n") != std::string::npos);
}

NOVA_TEST(pvd_collection_rejects_invalid_timesteps_without_partial_output) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const std::vector<std::vector<nova::io::VtkTimeSeriesEntry>> invalid{
        {},
        {{nan, "a.vtk"}},
        {{inf, "a.vtk"}},
        {{1.0, "a.vtk"}, {1.0, "b.vtk"}},
        {{2.0, "a.vtk"}, {1.0, "b.vtk"}}
    };
    for (const auto& entries : invalid) {
        std::ostringstream output;
        NOVA_REQUIRE(raises<std::invalid_argument>([&] {
            nova::io::write_pvd_collection(output, entries);
        }));
        NOVA_REQUIRE(output.str().empty());
    }
}

NOVA_TEST(pvd_collection_rejects_xml_injection_traversal_and_duplicates) {
    const std::vector<std::string> invalid{
        "../secret.vtk", "/tmp/a.vtk", "subdir/a.vtk",
        "evil&name.vtk", "bad\"quote.vtk", "<bad>.vtk",
        "bad\nname.vtk", ".hidden.vtk", "a..b.vtk",
        "a.vtu", "a.VTK", "a.vtk.pvd", ""
    };
    for (const auto& filename : invalid) {
        std::ostringstream output;
        NOVA_REQUIRE(raises<std::invalid_argument>([&] {
            nova::io::write_pvd_collection(output, {{0.0, filename}});
        }));
        NOVA_REQUIRE(output.str().empty());
    }
    std::ostringstream output;
    NOVA_REQUIRE(raises<std::invalid_argument>([&] {
        nova::io::write_pvd_collection(output, {{0.0, "same.vtk"}, {1.0, "same.vtk"}});
    }));
    NOVA_REQUIRE(output.str().empty());
}

NOVA_TEST(pvd_collection_rejects_nonwritable_stream) {
    std::ostringstream output;
    output.setstate(std::ios::badbit);
    NOVA_REQUIRE(raises<std::runtime_error>([&] {
        nova::io::write_pvd_collection(output, {{0.0, "step_0.vtk"}});
    }));
}
