#include "test_framework.hpp"
#include "nova/io/thermal_project.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace {
using namespace nova::model;

ThermalProject fixture() {
    ThermalProject p;
    p.geometry.lengths_m = {0.5, 0.25, 0.125};
    p.mesh.divisions = {6, 2, 1};
    p.material.conductivity_w_m_k = std::nextafter(2.0, 3.0);
    p.boundaries[BoxSide::XMin] = FixedTemperature{400};
    p.boundaries[BoxSide::XMax] = Convection{12, 300};
    p.solver = {1000, 1e-13, 1e-14, nova::linalg::PreconditionerType::Jacobi};
    return p;
}

class CommaDecimal : public std::numpunct<char> {
    char do_decimal_point() const override { return ','; }
};

template<class F> bool rejects(F operation) {
    try { operation(); }
    catch (const std::invalid_argument&) { return true; }
    return false;
}
} // namespace

NOVA_TEST(thermal_project_serialization_round_trips_in_si_independent_of_stream_locale) {
    auto p = fixture();
    p.boundaries[BoxSide::ZMax] = OutwardHeatFlux{-4.25};
    std::ostringstream out;
    out.imbue(std::locale(std::locale::classic(), new CommaDecimal));
    out << std::fixed << std::setprecision(1) << std::setw(12);
    nova::io::write_thermal_project(out, p);
    NOVA_REQUIRE(out.precision() == 1);
    NOVA_REQUIRE(out.width() == 12);
    std::istringstream input(out.str());
    input.imbue(out.getloc());
    const auto restored = nova::io::read_thermal_project(input);
    NOVA_REQUIRE(restored.material.conductivity_w_m_k == p.material.conductivity_w_m_k);
    NOVA_REQUIRE(restored.mesh.divisions == p.mesh.divisions);
    NOVA_REQUIRE(std::get<OutwardHeatFlux>(restored.boundaries.at(BoxSide::ZMax)).watts_per_m2 == -4.25);
    std::ostringstream second;
    nova::io::write_thermal_project(second, restored);
    NOVA_REQUIRE(second.str() == out.str());
}

NOVA_TEST(thermal_project_named_sides_keep_their_geometry_after_remeshing) {
    auto p = fixture();
    for (unsigned axis = 0; axis < 3; ++axis) {
        for (unsigned reverse = 0; reverse < 2; ++reverse) {
            p.boundaries.clear();
            p.boundaries[static_cast<BoxSide>(2 * axis + reverse)] = FixedTemperature{300};
            p.boundaries[static_cast<BoxSide>(2 * axis + 1 - reverse)] = OutwardHeatFlux{-20};
            for (const std::size_t divisions : {2U, 5U}) {
                p.mesh.divisions = {3, 2, 1};
                p.mesh.divisions[axis] = divisions;
                const auto result = solve_thermal_project(p);
                for (std::size_t node = 0; node < result.mesh.nodes.size(); ++node) {
                    const auto& point = result.mesh.nodes[node];
                    const double coordinate = std::array{point.x, point.y, point.z}[axis];
                    const double distance = reverse == 0 ? coordinate : p.geometry.lengths_m[axis] - coordinate;
                    NOVA_REQUIRE_NEAR(result.temperature.values[node], 300.0 + 10.0 * distance, 2e-8);
                }
                NOVA_REQUIRE_NEAR(result.balance.signed_balance_defect, 0, 2e-8);
            }
        }
    }
}

NOVA_TEST(thermal_project_validates_conflicting_references_physical_data_and_resource_limits) {
    auto p = fixture();
    p.material.id = p.geometry.id;
    NOVA_REQUIRE(rejects([&] { validate_thermal_project(p); }));
    p = fixture();
    p.geometry.id = "../box";
    NOVA_REQUIRE(rejects([&] { validate_thermal_project(p); }));
    p = fixture();
    p.mesh.divisions = {std::numeric_limits<std::size_t>::max(), 2, 4};
    NOVA_REQUIRE(rejects([&] { validate_thermal_project(p); }));
    p = fixture();
    p.boundaries[BoxSide::ZMin] = FixedTemperature{399};
    NOVA_REQUIRE(rejects([&] { validate_thermal_project(p); }));
    p.boundaries = {{BoxSide::XMin, OutwardHeatFlux{0}}};
    NOVA_REQUIRE(rejects([&] { validate_thermal_project(p); }));
    p.boundaries = {{BoxSide::ZMax, Convection{1, -273}}};
    NOVA_REQUIRE(rejects([&] { validate_thermal_project(p); }));
    p = fixture();
    p.solver.relative_tolerance = 1.0;
    NOVA_REQUIRE(rejects([&] { validate_thermal_project(p); }));
    p = fixture();
    p.generation_w_m3 = -1e8;
    NOVA_REQUIRE(rejects([&] { (void)solve_thermal_project(p); }));
}

NOVA_TEST(thermal_project_anchors_absolute_temperature_using_only_convection) {
    auto p = fixture();
    p.boundaries = {{BoxSide::YMax, Convection{10, 315}}};
    const auto result = solve_thermal_project(p);
    for (double value : result.temperature.values) NOVA_REQUIRE_NEAR(value, 315, 1e-8);
    NOVA_REQUIRE_NEAR(result.balance.robin_outward_heat, 0, 1e-9);
}
