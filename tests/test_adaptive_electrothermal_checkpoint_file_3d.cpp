#include "test_framework.hpp"

#include "nova/io/adaptive_electrothermal_checkpoint_file_3d.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace fs = std::filesystem;
using namespace nova::io;

class TestDirectory {
public:
    TestDirectory() {
        std::random_device entropy;
        for (unsigned i = 0; i < 64U; ++i) {
            std::ostringstream name;
            name << "nova-electrothermal-test-" << std::hex
                 << entropy() << '-' << entropy() << '-' << i;
            const auto candidate = fs::temp_directory_path() / name.str();
            if (fs::create_directory(candidate)) {
                directory_ = candidate;
                return;
            }
        }
        throw std::runtime_error("Could not reserve test directory");
    }
    ~TestDirectory() {
        std::error_code ignored;
        fs::remove_all(directory_, ignored);
    }
    [[nodiscard]] fs::path file() const { return directory_ / "restart.dat"; }
    [[nodiscard]] const fs::path& path() const { return directory_; }
private:
    fs::path directory_;
};

nova::fem::TetrahedralMesh3D unit_tetra() {
    return {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
            {{0, 1, 2, 3}}};
}

template<class F> bool rejects(F&& fn) {
    try { fn(); } catch (const std::exception&) { return true; }
    return false;
}

std::string bytes(const fs::path& path) {
    std::ifstream file(path, std::ios::binary);
    NOVA_REQUIRE(file.good());
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

void require_only_checkpoint_file(const TestDirectory& temp) {
    std::size_t count = 0;
    for (const auto& entry : fs::directory_iterator(temp.path())) {
        NOVA_REQUIRE(entry.path() == temp.file());
        ++count;
    }
    NOVA_REQUIRE(count == 1U);
}
} // namespace

NOVA_TEST(adaptive_electrothermal_checkpoint_file_3d_creates_and_replaces_losslessly) {
    TestDirectory temp;
    const auto mesh = unit_tetra();
    nova::physics::AdaptiveTransientElectrothermal3DOptions options;
    options.end_time = 1.0;
    const AdaptiveElectrothermalCheckpoint3D initial{
        std::nextafter(0.125, 1.0),
        {300.0, std::nextafter(301.0, 302.0), 302.5, 303.25},
        std::nextafter(0.05, 0.1)
    };
    save_adaptive_electrothermal_checkpoint_file_3d(
        temp.file(), mesh, "test_file_v1", options, initial);
    NOVA_REQUIRE(fs::is_regular_file(temp.file()));
    auto loaded = load_adaptive_electrothermal_checkpoint_file_3d(
        temp.file(), mesh, "test_file_v1", options);
    NOVA_REQUIRE(loaded.time == initial.time);
    NOVA_REQUIRE(loaded.temperature == initial.temperature);
    NOVA_REQUIRE(loaded.next_time_step == initial.next_time_step);
    const auto old_bytes = bytes(temp.file());

    const AdaptiveElectrothermalCheckpoint3D replacement{
        0.3, {305.0, 305.5, 306.25, 307.0}, 0.11
    };
    save_adaptive_electrothermal_checkpoint_file_3d(
        temp.file(), mesh, "test_file_v1", options, replacement);
    loaded = load_adaptive_electrothermal_checkpoint_file_3d(
        temp.file(), mesh, "test_file_v1", options);
    NOVA_REQUIRE(loaded.time == replacement.time);
    NOVA_REQUIRE(loaded.temperature == replacement.temperature);
    NOVA_REQUIRE(loaded.next_time_step == replacement.next_time_step);
    NOVA_REQUIRE(bytes(temp.file()) != old_bytes);
    require_only_checkpoint_file(temp);
}

NOVA_TEST(adaptive_electrothermal_checkpoint_file_3d_preserves_previous_on_validation_errors) {
    TestDirectory temp;
    const auto mesh = unit_tetra();
    nova::physics::AdaptiveTransientElectrothermal3DOptions options;
    options.end_time = 1.0;
    const AdaptiveElectrothermalCheckpoint3D initial{
        0.2, {300, 301, 302, 303}, 0.1
    };
    save_adaptive_electrothermal_checkpoint_file_3d(
        temp.file(), mesh, "test_file_v1", options, initial);
    const auto original = bytes(temp.file());

    auto invalid_state = initial;
    invalid_state.temperature[1] = std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects([&] {
        save_adaptive_electrothermal_checkpoint_file_3d(
            temp.file(), mesh, "test_file_v1", options, invalid_state);
    }));
    invalid_state = initial;
    invalid_state.next_time_step = -1.0;
    NOVA_REQUIRE(rejects([&] {
        save_adaptive_electrothermal_checkpoint_file_3d(
            temp.file(), mesh, "test_file_v1", options, invalid_state);
    }));
    auto invalid_options = options;
    invalid_options.maximum_attempts = 0U;
    NOVA_REQUIRE(rejects([&] {
        save_adaptive_electrothermal_checkpoint_file_3d(
            temp.file(), mesh, "test_file_v1", invalid_options, initial);
    }));
    NOVA_REQUIRE(rejects([&] {
        save_adaptive_electrothermal_checkpoint_file_3d(
            temp.file(), mesh, "invalid model id", options, initial);
    }));
    NOVA_REQUIRE(bytes(temp.file()) == original);
    require_only_checkpoint_file(temp);
}

NOVA_TEST(adaptive_electrothermal_checkpoint_file_3d_rejects_invalid_reads_and_targets) {
    TestDirectory temp;
    const auto mesh = unit_tetra();
    nova::physics::AdaptiveTransientElectrothermal3DOptions options;
    options.end_time = 1.0;
    const AdaptiveElectrothermalCheckpoint3D initial{
        0.2, {300, 301, 302, 303}, 0.1
    };
    NOVA_REQUIRE(rejects([&] {
        (void)load_adaptive_electrothermal_checkpoint_file_3d(
            temp.file(), mesh, "test_file_v1", options);
    }));
    save_adaptive_electrothermal_checkpoint_file_3d(
        temp.file(), mesh, "test_file_v1", options, initial);
    const auto original = bytes(temp.file());

    auto other_mesh = mesh;
    other_mesh.nodes[1].x = 2.0;
    NOVA_REQUIRE(rejects([&] {
        (void)load_adaptive_electrothermal_checkpoint_file_3d(
            temp.file(), other_mesh, "test_file_v1", options);
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)load_adaptive_electrothermal_checkpoint_file_3d(
            temp.file(), mesh, "different_model", options);
    }));
    NOVA_REQUIRE(rejects([&] {
        save_adaptive_electrothermal_checkpoint_file_3d(
            temp.path(), mesh, "test_file_v1", options, initial);
    }));
    const auto link = temp.path() / "checkpoint_link";
    std::error_code symlink_error;
    fs::create_symlink(temp.file(), link, symlink_error);
    if (!symlink_error) {
        NOVA_REQUIRE(rejects([&] {
            save_adaptive_electrothermal_checkpoint_file_3d(
                link, mesh, "test_file_v1", options, initial);
        }));
    }
    NOVA_REQUIRE(bytes(temp.file()) == original);

    {
        std::ofstream corrupt(temp.file(), std::ios::binary | std::ios::trunc);
        corrupt << "truncated checkpoint";
    }
    NOVA_REQUIRE(rejects([&] {
        (void)load_adaptive_electrothermal_checkpoint_file_3d(
            temp.file(), mesh, "test_file_v1", options);
    }));
}
