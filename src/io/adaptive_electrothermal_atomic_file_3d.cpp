#include "nova/io/adaptive_electrothermal_checkpoint_3d.hpp"

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>

namespace nova::io {
namespace {

// Directory creation is exclusive, avoiding exists-then-truncate races.
// A different staging directory is used for every writer. All staged files
// live on the destination filesystem to permit a single rename on publish.
class CheckpointStagingDirectory {
public:
    explicit CheckpointStagingDirectory(const std::filesystem::path& destination) {
        if (destination.empty() || destination.filename().empty()) {
            throw std::invalid_argument("Adaptive checkpoint destination must name a file.");
        }
        const auto parent = destination.has_parent_path()
            ? destination.parent_path() : std::filesystem::path(".");
        static std::atomic<unsigned long long> counter{0ULL};
        for (int i = 0; i < 64; ++i) {
            const auto tick =
                std::chrono::steady_clock::now().time_since_epoch().count();
            const auto candidate = parent /
                (destination.filename().string() + ".nova-stage-" +
                 std::to_string(tick) + "-" +
                 std::to_string(counter.fetch_add(1ULL)));
            std::error_code ec;
            if (std::filesystem::create_directory(candidate, ec)) {
                directory_ = candidate;
                return;
            }
            if (ec && ec != std::errc::file_exists) {
                throw std::filesystem::filesystem_error(
                    "Adaptive checkpoint staging directory creation failed",
                    candidate, ec);
            }
        }
        throw std::runtime_error("Could not reserve unique checkpoint staging directory.");
    }

    CheckpointStagingDirectory(const CheckpointStagingDirectory&) = delete;
    CheckpointStagingDirectory& operator=(const CheckpointStagingDirectory&) = delete;

    ~CheckpointStagingDirectory() {
        if (!directory_.empty()) {
            std::error_code ec;
            std::filesystem::remove_all(directory_, ec);
        }
    }

    std::filesystem::path staged_file() const {
        return directory_ / "checkpoint";
    }

private:
    std::filesystem::path directory_;
};

} // namespace

void write_adaptive_electrothermal_checkpoint_3d_file_atomic(
    const std::filesystem::path& destination,
    const fem::TetrahedralMesh3D& mesh,
    const std::string& model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& options,
    const AdaptiveElectrothermalCheckpoint3D& checkpoint
) {
    // Validate and serialize completely before opening or replacing files.
    std::ostringstream record;
    write_adaptive_electrothermal_checkpoint_3d(
        record, mesh, model_id, options, checkpoint);
    const std::string bytes = record.str();

    CheckpointStagingDirectory staging(destination);
    {
        std::ofstream out(staging.staged_file(),
                          std::ios::binary | std::ios::trunc);
        if (!out) {
            throw std::runtime_error("Cannot open staged adaptive checkpoint.");
        }
        out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
        out.flush();
        if (!out) {
            throw std::runtime_error("Cannot write staged adaptive checkpoint.");
        }
        out.close();
        if (!out) {
            throw std::runtime_error("Cannot close staged adaptive checkpoint.");
        }
    }
    // Replaces a regular file on filesystems that implement rename-over-file.
    // Failure leaves the prior destination intact. This is not fsync-durable.
    std::filesystem::rename(staging.staged_file(), destination);
}

AdaptiveElectrothermalCheckpoint3D
read_adaptive_electrothermal_checkpoint_3d_file(
    const std::filesystem::path& source,
    const fem::TetrahedralMesh3D& expected_mesh,
    const std::string& expected_model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& expected_options
) {
    std::ifstream input(source, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open adaptive electrothermal checkpoint.");
    }
    return read_adaptive_electrothermal_checkpoint_3d(
        input, expected_mesh, expected_model_id, expected_options);
}

} // namespace nova::io
