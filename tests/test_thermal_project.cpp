#include "test_framework.hpp"
#include "nova/io/thermal_project.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <string>

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



NOVA_TEST(thermal_project_schema_two_accepts_isotropic_material_and_canonicalizes_to_schema_one) {
    const auto project = fixture();
    std::ostringstream original;
    nova::io::write_thermal_project(original, project);
    auto schema_two = original.str();
    NOVA_REQUIRE(schema_two.rfind("NOVA_PROJECT 1\n", 0) == 0);
    schema_two.replace(0, 14, "NOVA_PROJECT 2");

    std::istringstream input(schema_two);
    const auto restored = nova::io::read_thermal_project(input);
    NOVA_REQUIRE(!restored.material.conductivity_tensor_w_m_k.has_value());
    NOVA_REQUIRE(restored.material.conductivity_w_m_k == project.material.conductivity_w_m_k);

    std::ostringstream canonical;
    nova::io::write_thermal_project(canonical, restored);
    NOVA_REQUIRE(canonical.str() == original.str());
}

NOVA_TEST(thermal_project_schema_two_round_trips_anisotropic_tensor_and_solves_exact_off_diagonal_case) {
    auto p = fixture();
    const nova::fem::SymmetricTensor3D full_tensor{2.0, 3.0, 4.0, 0.25, -0.15, 0.10};
    p.material.conductivity_tensor_w_m_k = full_tensor;

    std::ostringstream out;
    nova::io::write_thermal_project(out, p);
    NOVA_REQUIRE(out.str().find("NOVA_PROJECT 2\n") == 0);

    std::istringstream input(out.str());
    const auto restored = nova::io::read_thermal_project(input);
    NOVA_REQUIRE(restored.material.conductivity_tensor_w_m_k.has_value());
    const auto& round_trip = *restored.material.conductivity_tensor_w_m_k;
    NOVA_REQUIRE(round_trip.xx == full_tensor.xx);
    NOVA_REQUIRE(round_trip.yy == full_tensor.yy);
    NOVA_REQUIRE(round_trip.zz == full_tensor.zz);
    NOVA_REQUIRE(round_trip.xy == full_tensor.xy);
    NOVA_REQUIRE(round_trip.xz == full_tensor.xz);
    NOVA_REQUIRE(round_trip.yz == full_tensor.yz);
    std::ostringstream second;
    nova::io::write_thermal_project(second, restored);
    NOVA_REQUIRE(second.str() == out.str());

    // The full off-diagonal tensor has an exact linear field with prescribed
    // transverse Neumann fluxes and an analytical convective exit.
    auto anisotropic = restored;
    anisotropic.boundaries[BoxSide::YMin] = OutwardHeatFlux{-37.5};
    anisotropic.boundaries[BoxSide::YMax] = OutwardHeatFlux{37.5};
    anisotropic.boundaries[BoxSide::ZMin] = OutwardHeatFlux{22.5};
    anisotropic.boundaries[BoxSide::ZMax] = OutwardHeatFlux{-22.5};
    const auto result = solve_thermal_project(anisotropic);
    for (std::size_t node = 0; node < result.mesh.nodes.size(); ++node) {
        NOVA_REQUIRE_NEAR(result.temperature.values[node], 400.0 - 150.0 * result.mesh.nodes[node].x, 2e-8);
    }
    for (const auto& flux : result.cell_fluxes) {
        NOVA_REQUIRE_NEAR(flux.flux.x, 300.0, 2e-8);
        NOVA_REQUIRE_NEAR(flux.flux.y, 37.5, 2e-8);
        NOVA_REQUIRE_NEAR(flux.flux.z, -22.5, 2e-8);
    }
    NOVA_REQUIRE_NEAR(result.balance.robin_outward_heat, 9.375, 2e-8);
    NOVA_REQUIRE_NEAR(result.balance.neumann_outward_heat, 0.0, 2e-8);
    NOVA_REQUIRE_NEAR(result.balance.signed_balance_defect, 0.0, 2e-8);
}

NOVA_TEST(thermal_project_rejects_non_spd_anisotropic_tensors_and_schema_one_tensor_records) {
    auto p = fixture();
    for (const auto tensor : {
            nova::fem::SymmetricTensor3D{1.0, 1.0, 1.0, 2.0, 0.0, 0.0},
            nova::fem::SymmetricTensor3D{1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
            nova::fem::SymmetricTensor3D{1.0, 1.0, 1.0,
                std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}}) {
        p.material.conductivity_tensor_w_m_k = tensor;
        NOVA_REQUIRE(rejects([&] { validate_thermal_project(p); }));
    }

    const std::string legacy_tensor =
        "NOVA_PROJECT 1\n"
        "project p\n"
        "geometry g box 1 1 1\n"
        "mesh m g structured_tet 1 1 1\n"
        "material s anisotropic 2 3 4 0 0 0\n"
        "study st steady_heat m s\n"
        "source st 0\n"
        "solver st cg jacobi 1000 1e-10 1e-14\n"
        "boundary st g xmin temperature 300\n"
        "end\n";
    std::istringstream legacy_input(legacy_tensor);
    NOVA_REQUIRE(rejects([&] { (void)nova::io::read_thermal_project(legacy_input); }));
}
