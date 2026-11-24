#include "nova/io/adaptive_electrothermal_checkpoint_file_3d.hpp"

#include "atomic_file.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace nova::io {
namespace {
namespace fs = std::filesystem;

void require_file_path(const fs::path& path) {
    if (path.empty() || path.filename().empty()) {
        throw std::invalid_argument("Electrothermal checkpoint requires a filename.");
    }
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

    detail::write_atomic_file(destination, [&](std::ostream& output) {
        output << staged.view();
    }, [&](const fs::path& record) {
        const auto verified = load_adaptive_electrothermal_checkpoint_file_3d(
            record, mesh, model_id, options);
        if (verified.time != state.time || verified.next_time_step != state.next_time_step ||
            verified.temperature != state.temperature) {
            throw std::runtime_error("Electrothermal checkpoint read-back mismatch.");
        }
    });
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
