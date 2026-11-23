#include "test_framework.hpp"

#include "nova/io/adaptive_electrothermal_checkpoint_3d.hpp"
#include "nova/physics/adaptive_transient_electrothermal_3d.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <ios>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using namespace nova::fem;
using namespace nova::io;
using namespace nova::physics;

class CheckpointTestDirectory {
public:
    CheckpointTestDirectory() {
        const auto root = std::filesystem::temp_directory_path();
        for (int i = 0; i < 64; ++i) {
            auto p = root / ("nova-checkpoint-test-" +
                std::to_string(std::chrono::steady_clock::now()
                    .time_since_epoch().count()) + "-" + std::to_string(i));
            if (std::filesystem::create_directory(p)) {
                path_ = p;
                return;
            }
        }
        throw std::runtime_error("Cannot create isolated checkpoint test directory.");
    }
    ~CheckpointTestDirectory() {
        std::error_code ignored;
        if (!path_.empty()) std::filesystem::remove_all(path_, ignored);
    }
    const std::filesystem::path& path() const { return path_; }
private:
    std::filesystem::path path_;
};

class CommaDecimal final : public std::numpunct<char> {
protected:
    char do_decimal_point() const override { return ','; }
};

TetrahedralMesh3D unit_tetra() {
    return {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
            {{0, 1, 2, 3}}};
}

template<class F> bool rejects(F&& f) {
    try { f(); } catch (const std::exception&) { return true; }
    return false;
}

std::string write_test_record(
    const TetrahedralMesh3D& mesh,
    const AdaptiveTransientElectrothermal3DOptions& options,
    const AdaptiveElectrothermalCheckpoint3D& state
) {
    std::ostringstream out;
    write_adaptive_electrothermal_checkpoint_3d(
        out, mesh, "unit_cube_v1", options, state);
    return out.str();
}

std::string replaced(std::string text,
                     const std::string& old_value,
                     const std::string& new_value) {
    const auto offset = text.find(old_value);
    NOVA_REQUIRE(offset != std::string::npos);
    text.replace(offset, old_value.size(), new_value);
    return text;
}

std::string change_field(std::string text, const std::string& label,
                         const std::string& value) {
    const auto first = text.find(label);
    NOVA_REQUIRE(first != std::string::npos);
    const auto end = text.find('\n', first);
    NOVA_REQUIRE(end != std::string::npos);
    text.replace(first, end - first, label + value);
    return text;
}

bool face_at_x(const TetrahedralMesh3D& mesh,
               const BoundaryFace3D& face, double x) {
    return std::all_of(face.begin(), face.end(),
        [&](std::size_t i) { return mesh.nodes[i].x == x; });
}

TransientElectrothermal3DProblem make_uniform_ode() {
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 1, 1, 1);
    TransientElectrothermal3DProblem p;
    p.electrical.mesh = box.mesh;
    p.electrical.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    p.thermal.diffusion.mesh = box.mesh;
    p.thermal.diffusion.cell_diffusivity.assign(box.mesh.elements.size(), 4.0);
    p.thermal.cell_volumetric_capacity.assign(box.mesh.elements.size(), 2.0);
    p.thermal.initial_temperature.assign(box.mesh.nodes.size(), 300.0);
    for (const auto& face : box.boundary_faces) {
        if (face_at_x(box.mesh, face, 0.0) ||
            face_at_x(box.mesh, face, 1.0)) {
            p.electrical.dirichlet_faces.push_back({
                face, [](double x, double, double) { return 10.0 * x; }
            });
        }
    }
    p.electrical_conductivity_scale = [](double t) {
        return 1.0 + 0.01 * (t - 300.0);
    };
    p.thermal_conductivity_scale = [](double t) {
        return 1.0 + 0.01 * (t - 300.0);
    };
    return p;
}
} // namespace

NOVA_TEST(adaptive_electrothermal_checkpoint_3d_roundtrip_exact_and_stream_safe) {
    const auto mesh = unit_tetra();
    AdaptiveTransientElectrothermal3DOptions options;
    options.end_time = 1.0;
    options.coupling.thermal_time.solver_options.relative_tolerance = 1e-12;
    options.coupling.electrical_solver_options.absolute_tolerance = 1e-15;
    const AdaptiveElectrothermalCheckpoint3D state{
        std::nextafter(0.125, 1.0),
        {300.5, std::nextafter(1.0, 2.0), -1e-200, 1.2345678901234567e100},
        std::nextafter(0.1, 0.2)
    };
    std::ostringstream output;
    output.imbue(std::locale(std::locale::classic(), new CommaDecimal));
    output << std::fixed << std::setprecision(2);
    const auto flags = output.flags();
    write_adaptive_electrothermal_checkpoint_3d(
        output, mesh, "unit_cube_v1", options, state);
    NOVA_REQUIRE(output.flags() == flags);
    NOVA_REQUIRE(output.precision() == 2);
    NOVA_REQUIRE(output.str().find("version 1") != std::string::npos);
    NOVA_REQUIRE(output.str().find("thermal_state\nNOVA_THERMAL_CHECKPOINT_3D 1") !=
                 std::string::npos);

    std::istringstream input(output.str());
    input.imbue(std::locale(std::locale::classic(), new CommaDecimal));
    const auto loaded = read_adaptive_electrothermal_checkpoint_3d(
        input, mesh, "unit_cube_v1", options);
    NOVA_REQUIRE(loaded.time == state.time);
    NOVA_REQUIRE(loaded.next_time_step == state.next_time_step);
    NOVA_REQUIRE(loaded.temperature == state.temperature);
}

NOVA_TEST(adaptive_electrothermal_checkpoint_3d_rejects_mismatched_options_model_and_mesh) {
    const auto mesh = unit_tetra();
    AdaptiveTransientElectrothermal3DOptions options;
    options.end_time = 0.5;
    const AdaptiveElectrothermalCheckpoint3D state{
        0.1, {300, 301, 302, 303}, 0.05
    };
    const auto good = write_test_record(mesh, options, state);
    auto different = options;
    different.relative_error *= 2.0;
    std::istringstream changed_tolerance(good);
    NOVA_REQUIRE(rejects([&] {
        (void)read_adaptive_electrothermal_checkpoint_3d(
            changed_tolerance, mesh, "unit_cube_v1", different);
    }));
    different = options;
    different.coupling.relaxation = 0.5;
    std::istringstream changed_coupling(good);
    NOVA_REQUIRE(rejects([&] {
        (void)read_adaptive_electrothermal_checkpoint_3d(
            changed_coupling, mesh, "unit_cube_v1", different);
    }));
    different = options;
    different.coupling.thermal_time.solver_options.max_iterations = 777;
    std::istringstream changed_solver(good);
    NOVA_REQUIRE(rejects([&] {
        (void)read_adaptive_electrothermal_checkpoint_3d(
            changed_solver, mesh, "unit_cube_v1", different);
    }));
    std::istringstream changed_id(good);
    NOVA_REQUIRE(rejects([&] {
        (void)read_adaptive_electrothermal_checkpoint_3d(
            changed_id, mesh, "some_other_model", options);
    }));
    auto altered_mesh = mesh;
    altered_mesh.nodes[1].x += 0.1;
    std::istringstream changed_mesh(good);
    NOVA_REQUIRE(rejects([&] {
        (void)read_adaptive_electrothermal_checkpoint_3d(
            changed_mesh, altered_mesh, "unit_cube_v1", options);
    }));
    altered_mesh = mesh;
    std::swap(altered_mesh.elements[0][1], altered_mesh.elements[0][2]);
    std::istringstream changed_connectivity(good);
    NOVA_REQUIRE(rejects([&] {
        (void)read_adaptive_electrothermal_checkpoint_3d(
            changed_connectivity, altered_mesh, "unit_cube_v1", options);
    }));
}

NOVA_TEST(adaptive_electrothermal_checkpoint_3d_rejects_bad_records_without_partial_write) {
    const auto mesh = unit_tetra();
    AdaptiveTransientElectrothermal3DOptions options;
    options.end_time = 0.5;
    const AdaptiveElectrothermalCheckpoint3D state{
        0.1, {300, 301, 302, 303}, 0.05
    };
    const std::string good = write_test_record(mesh, options, state);
    const std::vector<std::string> corrupt{
        "",
        replaced(good, "version 1", "version 2"),
        change_field(good, "min_dt ", "-1"),
        change_field(good, "next_dt ", "nan"),
        replaced(good, "thermal_state\n", "thermal_state\nBROKEN\n"),
        good.substr(0, good.size() - 4),
        good + "trailing\n"
    };
    for (const auto& record : corrupt) {
        std::istringstream in(record);
        NOVA_REQUIRE(rejects([&] {
            (void)read_adaptive_electrothermal_checkpoint_3d(
                in, mesh, "unit_cube_v1", options);
        }));
    }
    auto bad = state;
    bad.next_time_step = std::numeric_limits<double>::infinity();
    std::ostringstream no_partial;
    NOVA_REQUIRE(rejects([&] {
        write_adaptive_electrothermal_checkpoint_3d(
            no_partial, mesh, "unit_cube_v1", options, bad);
    }));
    NOVA_REQUIRE(no_partial.str().empty());
    bad = state;
    bad.temperature[0] = std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects([&] {
        write_adaptive_electrothermal_checkpoint_3d(
            no_partial, mesh, "unit_cube_v1", options, bad);
    }));
    NOVA_REQUIRE(no_partial.str().empty());
    NOVA_REQUIRE(rejects([&] {
        write_adaptive_electrothermal_checkpoint_3d(
            no_partial, mesh, "invalid model id", options, state);
    }));
    NOVA_REQUIRE(no_partial.str().empty());
}

NOVA_TEST(adaptive_electrothermal_checkpoint_3d_nonlinear_restart_matches_continuous_run) {
    const auto problem = make_uniform_ode();
    AdaptiveTransientElectrothermal3DOptions options;
    options.end_time = 0.3;
    options.initial_time_step = 0.3;
    options.maximum_time_step = 0.3;
    options.minimum_time_step = 1e-5;
    options.absolute_error_kelvin = 0.05;
    options.relative_error = 0.0;
    std::ostringstream output;
    std::size_t callbacks = 0U;
    options.on_accepted_step = [&](const AdaptiveTransientElectrothermal3DStep& step,
                                   const std::vector<double>& field) {
        ++callbacks;
        if (callbacks == 1U) {
            write_adaptive_electrothermal_checkpoint_3d(
                output, problem.electrical.mesh, "unit_cube_v1", options,
                {step.time, field, step.next_time_step});
        }
    };
    const auto continuous = solve_adaptive_transient_electrothermal_3d(problem, options);
    NOVA_REQUIRE(continuous.completed);
    NOVA_REQUIRE(continuous.rejected_steps >= 1U);
    NOVA_REQUIRE(callbacks == continuous.accepted_steps.size());
    NOVA_REQUIRE(callbacks >= 3U);

    std::istringstream input(output.str());
    const auto state = read_adaptive_electrothermal_checkpoint_3d(
        input, problem.electrical.mesh, "unit_cube_v1", options);
    NOVA_REQUIRE_NEAR(state.time, continuous.accepted_steps.front().time, 0.0);
    auto resumed_problem = problem;
    resumed_problem.thermal.initial_temperature = state.temperature;
    auto resumed_options = options;
    resumed_options.coupling.thermal_time.start_time = state.time;
    resumed_options.initial_time_step = state.next_time_step;
    resumed_options.on_accepted_step = {};
    const auto resumed = solve_adaptive_transient_electrothermal_3d(
        resumed_problem, resumed_options);
    NOVA_REQUIRE(resumed.completed);
    NOVA_REQUIRE_NEAR(resumed.final_time, continuous.final_time, 1e-14);
    NOVA_REQUIRE(resumed.accepted_steps.size() + 1U ==
                 continuous.accepted_steps.size());
    for (std::size_t i = 0; i < resumed.accepted_steps.size(); ++i) {
        NOVA_REQUIRE_NEAR(resumed.accepted_steps[i].time,
                          continuous.accepted_steps[i + 1U].time, 1e-12);
    }
    for (std::size_t i = 0; i < resumed.temperature.size(); ++i) {
        NOVA_REQUIRE_NEAR(resumed.temperature[i], continuous.temperature[i], 1e-8);
    }
}


NOVA_TEST(adaptive_electrothermal_checkpoint_3d_atomic_file_replacement) {
    CheckpointTestDirectory directory;
    const auto path = directory.path() / "restart.dat";
    const auto mesh = unit_tetra();
    AdaptiveTransientElectrothermal3DOptions options;
    options.end_time = 1.0;
    const AdaptiveElectrothermalCheckpoint3D before{0.2, {300, 301, 302, 303}, 0.05};
    const AdaptiveElectrothermalCheckpoint3D after{0.4, {305, 306, 307, 308}, 0.1};

    write_adaptive_electrothermal_checkpoint_3d_file_atomic(
        path, mesh, "unit_cube_v1", options, before);
    auto loaded = read_adaptive_electrothermal_checkpoint_3d_file(
        path, mesh, "unit_cube_v1", options);
    NOVA_REQUIRE(loaded.time == before.time);
    NOVA_REQUIRE(loaded.temperature == before.temperature);
    NOVA_REQUIRE(loaded.next_time_step == before.next_time_step);

    write_adaptive_electrothermal_checkpoint_3d_file_atomic(
        path, mesh, "unit_cube_v1", options, after);
    loaded = read_adaptive_electrothermal_checkpoint_3d_file(
        path, mesh, "unit_cube_v1", options);
    NOVA_REQUIRE(loaded.time == after.time);
    NOVA_REQUIRE(loaded.temperature == after.temperature);
    NOVA_REQUIRE(loaded.next_time_step == after.next_time_step);

    std::size_t entries = 0U;
    for (const auto& entry : std::filesystem::directory_iterator(directory.path())) {
        NOVA_REQUIRE(entry.path() == path);
        ++entries;
    }
    NOVA_REQUIRE(entries == 1U);
}

NOVA_TEST(adaptive_electrothermal_checkpoint_3d_failed_publish_keeps_old_file) {
    CheckpointTestDirectory directory;
    const auto path = directory.path() / "restart.dat";
    const auto mesh = unit_tetra();
    AdaptiveTransientElectrothermal3DOptions options;
    options.end_time = 1.0;
    const AdaptiveElectrothermalCheckpoint3D state{0.2, {300, 301, 302, 303}, 0.05};

    write_adaptive_electrothermal_checkpoint_3d_file_atomic(
        path, mesh, "unit_cube_v1", options, state);
    auto invalid = state;
    invalid.temperature[0] = std::numeric_limits<double>::infinity();
    NOVA_REQUIRE(rejects([&] {
        write_adaptive_electrothermal_checkpoint_3d_file_atomic(
            path, mesh, "unit_cube_v1", options, invalid);
    }));
    NOVA_REQUIRE(rejects([&] {
        write_adaptive_electrothermal_checkpoint_3d_file_atomic(
            path, mesh, "invalid model id", options, state);
    }));

    // Publishing a file over a nonempty directory is invalid on both POSIX
    // and Windows. Confirm cleanup and untouched unrelated contents.
    const auto blocked = directory.path() / "occupied";
    std::filesystem::create_directory(blocked);
    {
        std::ofstream sentinel(blocked / "keep.txt");
        sentinel << "untouched";
    }
    NOVA_REQUIRE(rejects([&] {
        write_adaptive_electrothermal_checkpoint_3d_file_atomic(
            blocked, mesh, "unit_cube_v1", options, state);
    }));
    NOVA_REQUIRE(std::filesystem::exists(blocked / "keep.txt"));
    const auto restored = read_adaptive_electrothermal_checkpoint_3d_file(
        path, mesh, "unit_cube_v1", options);
    NOVA_REQUIRE(restored.time == state.time);
    NOVA_REQUIRE(restored.temperature == state.temperature);
    NOVA_REQUIRE(restored.next_time_step == state.next_time_step);
    std::size_t entries = 0U;
    for (const auto& entry : std::filesystem::directory_iterator(directory.path())) {
        NOVA_REQUIRE(entry.path() == path || entry.path() == blocked);
        ++entries;
    }
    NOVA_REQUIRE(entries == 2U);
    NOVA_REQUIRE(rejects([&] {
        (void)read_adaptive_electrothermal_checkpoint_3d_file(
            directory.path() / "missing.dat", mesh, "unit_cube_v1", options);
    }));
}
