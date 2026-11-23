#pragma once

#include "nova/io/adaptive_electrothermal_checkpoint_3d.hpp"

#include <filesystem>
#include <string>

namespace nova::io {

// Serialize and validate the accepted electrothermal state before touching
// disk. Write and close a uniquely reserved staging file in the destination
// directory, verify its contents, then replace the destination.
// A failed validation/write/replacement preserves the previously committed
// file under normal filesystem error semantics. One writer per path is
// required; this is NOT a power-loss durable or cross-process lock protocol.
// POSIX same-filesystem rename provides atomic name replacement; Windows
// uses ReplaceFileW for an existing file, MoveFileExW for a new file.
// A parent directory must exist. Symlink destinations are rejected.
void save_adaptive_electrothermal_checkpoint_file_3d(
    const std::filesystem::path& destination,
    const fem::TetrahedralMesh3D& mesh,
    const std::string& model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& options,
    const AdaptiveElectrothermalCheckpoint3D& state
);

// Open and fully validate an on-disk checkpoint against the expected mesh,
// physical model token and controller/solver settings. Never updates the
// caller's state or model during a failed read.
[[nodiscard]] AdaptiveElectrothermalCheckpoint3D
load_adaptive_electrothermal_checkpoint_file_3d(
    const std::filesystem::path& source,
    const fem::TetrahedralMesh3D& expected_mesh,
    const std::string& expected_model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& expected_options
);

} // namespace nova::io
