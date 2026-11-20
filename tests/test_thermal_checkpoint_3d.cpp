#include "test_framework.hpp"

#include "nova/fem/transient_heat_3d.hpp"
#include "nova/io/thermal_checkpoint_3d.hpp"

#include <cmath>
#include <iomanip>
#include <ios>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

namespace {

class CommaDecimal final : public std::numpunct<char> {
protected:
    char do_decimal_point() const override { return ','; }
};

nova::fem::TetrahedralMesh3D unit_mesh() {
    return {
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
        {{0, 1, 2, 3}}
    };
}

template<class Exception, class F>
bool raises(F&& action) {
    try {
        action();
    } catch (const Exception&) {
        return true;
    }
    return false;
}

std::string checkpoint_text(
    const nova::fem::TetrahedralMesh3D& mesh,
    double time,
    const std::vector<double>& values
) {
    std::ostringstream output;
    nova::io::write_transient_heat_checkpoint_3d(output, mesh, time, values);
    return output.str();
}

} // namespace

NOVA_TEST(thermal_checkpoint_3d_roundtrips_exact_doubles_and_preserves_stream_format) {
    const auto mesh = unit_mesh();
    const std::vector<double> values{
        300.5, std::nextafter(1.0, 2.0), -1e-200, 1.2345678901234567e100
    };
    const std::locale comma_locale(std::locale::classic(), new CommaDecimal);
    const double time = std::nextafter(3.125, 4.0);
    std::ostringstream output;
    output.imbue(comma_locale);
    output << std::fixed << std::setprecision(2);
    const auto flags = output.flags();
    nova::io::write_transient_heat_checkpoint_3d(output, mesh, time, values);
    NOVA_REQUIRE(output.flags() == flags);
    NOVA_REQUIRE(output.precision() == 2);
    NOVA_REQUIRE(output.str().find("time 3.") != std::string::npos);

    std::istringstream input(output.str());
    input.imbue(comma_locale);
    const auto loaded = nova::io::read_transient_heat_checkpoint_3d(input, mesh);
    NOVA_REQUIRE(loaded.time == time);
    NOVA_REQUIRE(loaded.temperature == values);
}

NOVA_TEST(thermal_checkpoint_3d_rejects_nonfinite_data_before_writing) {
    const auto mesh = unit_mesh();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (double bad : {nan, inf, -inf}) {
        std::ostringstream bad_time;
        NOVA_REQUIRE(raises<std::invalid_argument>([&] {
            nova::io::write_transient_heat_checkpoint_3d(
                bad_time, mesh, bad, {300, 301, 302, 303});
        }));
        NOVA_REQUIRE(bad_time.str().empty());

        std::ostringstream bad_field;
        NOVA_REQUIRE(raises<std::invalid_argument>([&] {
            nova::io::write_transient_heat_checkpoint_3d(
                bad_field, mesh, 0.5, {300, bad, 302, 303});
        }));
        NOVA_REQUIRE(bad_field.str().empty());
    }

    std::ostringstream wrong_size;
    NOVA_REQUIRE(raises<std::invalid_argument>([&] {
        nova::io::write_transient_heat_checkpoint_3d(
            wrong_size, mesh, 0.5, {300, 301, 302});
    }));
    NOVA_REQUIRE(wrong_size.str().empty());

    std::ostringstream blocked;
    blocked.setstate(std::ios::badbit);
    NOVA_REQUIRE(raises<std::runtime_error>([&] {
        nova::io::write_transient_heat_checkpoint_3d(
            blocked, mesh, 0.5, {300, 301, 302, 303});
    }));
}

NOVA_TEST(thermal_checkpoint_3d_rejects_different_mesh_with_identical_counts) {
    const auto mesh = unit_mesh();
    const auto text = checkpoint_text(mesh, 0.5, {300, 301, 302, 303});
    auto modified = mesh;
    modified.nodes[1].x = 1.25;
    std::istringstream geometry_input(text);
    NOVA_REQUIRE(raises<std::invalid_argument>([&] {
        (void)nova::io::read_transient_heat_checkpoint_3d(geometry_input, modified);
    }));
    modified = mesh;
    std::swap(modified.elements[0][1], modified.elements[0][2]);
    std::istringstream connectivity_input(text);
    NOVA_REQUIRE(raises<std::invalid_argument>([&] {
        (void)nova::io::read_transient_heat_checkpoint_3d(connectivity_input, modified);
    }));
}

NOVA_TEST(thermal_checkpoint_3d_rejects_corruption_version_and_trailing_data) {
    const auto mesh = unit_mesh();
    const std::string good = checkpoint_text(mesh, 4.25, {300, 301, 302, 303});
    const std::vector<std::string> malformed{
        "",
        "NOVA_THERMAL_CHECKPOINT_3D 2\n" + good.substr(good.find('\n') + 1),
        good.substr(0, good.size() - 4),
        good + "unexpected\n",
        good.substr(0, good.find("nodes 4")) +
            "nodes 999999999\n" + good.substr(good.find("nodes 4") + 8),
        good.substr(0, good.find("300\n")) + "nan\n" +
            good.substr(good.find("300\n") + 4)
    };
    for (const auto& text : malformed) {
        std::istringstream input(text);
        NOVA_REQUIRE(raises<std::invalid_argument>([&] {
            (void)nova::io::read_transient_heat_checkpoint_3d(input, mesh);
        }));
    }
}

NOVA_TEST(thermal_checkpoint_3d_restores_transient_solver_without_nodal_drift) {
    using namespace nova::fem;
    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                              TransientHeat3DScheme::CrankNicolson}) {
        TransientHeat3DProblem problem;
        problem.diffusion.mesh = unit_mesh();
        problem.diffusion.cell_diffusivity = {1.0};
        problem.cell_volumetric_capacity = {2.0};
        problem.initial_temperature = {0.0, 0.0, 0.0, 10.0};
        problem.diffusion.dirichlet_constraints = {
            {0, 0.0}, {1, 0.0}, {2, 0.0}
        };
        TransientHeat3DOptions options;
        options.scheme = scheme;
        options.start_time = 1.25;
        options.time_step = 0.05;
        options.steps = 6;
        options.solver_options.relative_tolerance = 1e-12;

        const auto continuous = solve_transient_heat_3d(problem, options);
        options.steps = 2;
        const auto first = solve_transient_heat_3d(problem, options);
        std::istringstream saved(checkpoint_text(
            problem.diffusion.mesh, first.final_time, first.temperature));
        const auto loaded = nova::io::read_transient_heat_checkpoint_3d(
            saved, problem.diffusion.mesh);

        problem.initial_temperature = loaded.temperature;
        options.start_time = loaded.time;
        options.steps = 4;
        const auto resumed = solve_transient_heat_3d(problem, options);
        NOVA_REQUIRE_NEAR(resumed.final_time, continuous.final_time, 1e-12);
        NOVA_REQUIRE(resumed.temperature.size() == continuous.temperature.size());
        for (std::size_t i = 0; i < resumed.temperature.size(); ++i) {
            NOVA_REQUIRE_NEAR(
                resumed.temperature[i], continuous.temperature[i], 1e-11);
        }
    }
}
