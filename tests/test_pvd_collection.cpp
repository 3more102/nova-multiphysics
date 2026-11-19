#include "test_framework.hpp"
#include "nova/io/pvd_collection.hpp"

#include <ios>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

template<class Exception, class Action>
bool throws_as(Action&& action) {
    try { action(); } catch (const Exception&) { return true; }
    return false;
}

} // namespace

NOVA_TEST(pvd_collection_writes_ordered_time_frames_at_double_precision) {
    std::ostringstream output;
    output.precision(4);
    const double t = 0.12345678901234568;
    const std::vector<nova::io::ParaviewTimeFrame> frames{
        {0.0, "cooling_000000.vtk"},
        {t, "cooling_000001.vtk"},
        {0.5, "cooling_000002.vtk"}
    };
    nova::io::write_paraview_collection(output, frames);
    NOVA_REQUIRE(output.precision() == 4);
    const auto xml = output.str();
    NOVA_REQUIRE(xml.find("<?xml version=\"1.0\"?>\n") == 0U);
    NOVA_REQUIRE(xml.find("<VTKFile type=\"Collection\"") != std::string::npos);
    NOVA_REQUIRE(xml.find("timestep=\"0\" group=\"\" part=\"0\" file=\"cooling_000000.vtk\"")
                 != std::string::npos);
    NOVA_REQUIRE(xml.find("timestep=\"0.12345678901234568\"") != std::string::npos);
    NOVA_REQUIRE(xml.find("timestep=\"0.5\" group=\"\" part=\"0\" file=\"cooling_000002.vtk\"")
                 != std::string::npos);
    NOVA_REQUIRE(xml.find("</Collection>\n</VTKFile>\n") != std::string::npos);
}

NOVA_TEST(pvd_collection_rejects_nonmonotone_nan_and_unsafe_frames_before_writing) {
    using nova::io::ParaviewTimeFrame;
    std::ostringstream output;
    const auto rejects = [&](const std::vector<ParaviewTimeFrame>& frames) {
        NOVA_REQUIRE(throws_as<std::invalid_argument>([&] {
            nova::io::write_paraview_collection(output, frames);
        }));
        NOVA_REQUIRE(output.str().empty());
    };
    rejects({});
    rejects({{0., "a.vtk"}, {0., "b.vtk"}});
    rejects({{1., "a.vtk"}, {0., "b.vtk"}});
    rejects({{std::numeric_limits<double>::quiet_NaN(), "a.vtk"}});
    rejects({{std::numeric_limits<double>::infinity(), "a.vtk"}});
    rejects({{0., "../escape.vtk"}});
    rejects({{0., "bad\"name.vtk"}});
    rejects({{0., "bad&name.vtk"}});
    rejects({{0., "bad\\name.vtk"}});
    rejects({{0., "."}});
}

NOVA_TEST(pvd_collection_rejects_unwritable_stream) {
    std::ostringstream output;
    output.setstate(std::ios::badbit);
    NOVA_REQUIRE(throws_as<std::runtime_error>([&] {
        nova::io::write_paraview_collection(output, {{0., "snapshot.vtk"}});
    }));
}
