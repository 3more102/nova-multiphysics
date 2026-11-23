#include "nova/io/adaptive_electrothermal_checkpoint_file_3d.hpp"

#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace nova::io {
namespace {
namespace fs = std::filesystem;

void require_file_path(const fs::path& path) {
    if (path.empty() || path.filename().empty()) {
        throw std::invalid_argument("Electrothermal checkpoint requires a filename.");
    }
}

void reject_unsafe_destination(const fs::path& path) {
    const auto status = fs::symlink_status(path);
    if (fs::is_symlink(status) ||
        (fs::exists(status) && !fs::is_regular_file(status))) {
        throw std::invalid_argument(
            "Electrothermal checkpoint destination must not be a symlink or non-file.");
    }
}

fs::path reserve_staging_directory(const fs::path& parent) {
    std::random_device randomness;
    for (unsigned attempt = 0; attempt < 64U; ++attempt) {
        std::ostringstream name;
        name << ".nova-electrothermal-" << std::hex
             << randomness() << '-' << randomness() << '-' << attempt;
        const auto candidate = parent / name.str();
        if (fs::create_directory(candidate)) {
            return candidate;
        }
    }
    throw std::runtime_error("Cannot reserve a unique electrothermal checkpoint staging directory.");
}

class StagingDirectory {
public:
    explicit StagingDirectory(fs::path directory)
        : directory_(std::move(directory)) {}

    StagingDirectory(const StagingDirectory&) = delete;
    StagingDirectory& operator=(const StagingDirectory&) = delete;

    ~StagingDirectory() {
        std::error_code ignored;
        fs::remove_all(directory_, ignored);
    }

    [[nodiscard]] fs::path record() const {
        return directory_ / "checkpoint.dat";
    }

private:
    fs::path directory_;
};

void replace_file(const fs::path& staged, const fs::path& destination) {
    reject_unsafe_destination(destination);
#ifdef _WIN32
    // ReplaceFileW replaces an existing file without an application-level
    // remove gap. MoveFileExW handles first-time creation. Both files share
    // a parent filesystem; neither call promises power-failure durability.
    const bool existing = fs::exists(destination);
    const BOOL success = existing
        ? ::ReplaceFileW(destination.c_str(), staged.c_str(), nullptr,
                         0, nullptr, nullptr)
        : ::MoveFileExW(staged.c_str(), destination.c_str(),
                        MOVEFILE_WRITE_THROUGH);
    if (!success) {
        throw std::system_error(
            static_cast<int>(::GetLastError()), std::system_category(),
            "Failed to replace electrothermal checkpoint");
    }
#else
    fs::rename(staged, destination);
#endif
}
} // namespace

void save_adaptive_electrothermal_checkpoint_file_3d(
    const std::filesystem::path& destination,
    const fem::TetrahedralMesh3D& mesh,
    const std::string& model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& options,
    const AdaptiveElectrothermalCheckpoint3D& state
) {
    require_file_path(destination);
    // All model, mesh, settings and temperature validation happens before
    // creating a temporary file or changing the previous checkpoint.
    std::ostringstream staged;
    write_adaptive_electrothermal_checkpoint_3d(
        staged, mesh, model_id, options, state);

    reject_unsafe_destination(destination);
    const fs::path parent = destination.parent_path().empty()
        ? fs::path(".") : destination.parent_path();
    StagingDirectory temporary(reserve_staging_directory(parent));
    const fs::path record_path = temporary.record();
    {
        std::ofstream file(record_path, std::ios::binary | std::ios::trunc);
        if (!file) {
            throw std::runtime_error("Cannot open electrothermal checkpoint staging file.");
        }
        file << staged.str();
        file.flush();
        if (!file) {
            throw std::runtime_error("Cannot flush electrothermal checkpoint staging file.");
        }
        file.close();
        if (!file) {
            throw std::runtime_error("Cannot close electrothermal checkpoint staging file.");
        }
    }
    // Reject an incomplete/damaged staging record before any replacement.
    const auto verified = load_adaptive_electrothermal_checkpoint_file_3d(
        record_path, mesh, model_id, options);
    if (verified.time != state.time ||
        verified.next_time_step != state.next_time_step ||
        verified.temperature != state.temperature) {
        throw std::runtime_error("Electrothermal checkpoint read-back mismatch.");
    }
    replace_file(record_path, destination);
}

AdaptiveElectrothermalCheckpoint3D
load_adaptive_electrothermal_checkpoint_file_3d(
    const std::filesystem::path& source,
    const fem::TetrahedralMesh3D& expected_mesh,
    const std::string& expected_model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& expected_options
) {
    require_file_path(source);
    std::ifstream file(source, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open electrothermal checkpoint file.");
    }
    return read_adaptive_electrothermal_checkpoint_3d(
        file, expected_mesh, expected_model_id, expected_options);
}

} // namespace nova::io
