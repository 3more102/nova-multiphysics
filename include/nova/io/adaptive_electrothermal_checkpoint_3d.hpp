#pragma once

#include "nova/fem/mesh3d.hpp"
#include "nova/physics/adaptive_transient_electrothermal_3d.hpp"

#include <filesystem>
#include <iosfwd>
#include <string>
#include <vector>

namespace nova::io {

// Last committed adaptive two-half-step state only. The proposed timestep
// is retained so a restarted run uses the same controller proposal.
struct AdaptiveElectrothermalCheckpoint3D {
    double time{};
    std::vector<double> temperature;
    double next_time_step{};
};

// Single versioned text record: validated mesh and temperatures, physical
// time, next timestep, caller-assigned model identity and all numerical
// settings that affect the adaptive controller and its nested solvers.
// Input model functions, source/BC callbacks and full attempt history are
// NOT serializable. model_id is an external identity token, not a computed
// physics checksum. The caller must recreate an equivalent physical model.
//
// The writer stages data before output, validates numerical state and does
// not change the caller's stream formatting. Application code is responsible
// for durable/atomic file replacement and checking stream closure.
void write_adaptive_electrothermal_checkpoint_3d(
    std::ostream& output,
    const fem::TetrahedralMesh3D& mesh,
    const std::string& model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& options,
    const AdaptiveElectrothermalCheckpoint3D& checkpoint
);

// Exact comparison against the caller's supplied numerical options,
// excluding the moving start time, initial dt and overwritten single-step
// working fields. Requires the same caller-provided model_id and mesh.
// Rejects format/version mismatch, truncated/extra data and nonfinite state.
[[nodiscard]] AdaptiveElectrothermalCheckpoint3D
read_adaptive_electrothermal_checkpoint_3d(
    std::istream& input,
    const fem::TetrahedralMesh3D& expected_mesh,
    const std::string& expected_model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& expected_options
);

// Path API: publish a fully written sibling staged file through rename.
// Replacement atomicity depends on the platform/filesystem. No fsync or
// cross-process writer arbitration is performed.
void write_adaptive_electrothermal_checkpoint_3d_file_atomic(
    const std::filesystem::path& destination,
    const fem::TetrahedralMesh3D& mesh,
    const std::string& model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& options,
    const AdaptiveElectrothermalCheckpoint3D& checkpoint
);

[[nodiscard]] AdaptiveElectrothermalCheckpoint3D
read_adaptive_electrothermal_checkpoint_3d_file(
    const std::filesystem::path& source,
    const fem::TetrahedralMesh3D& expected_mesh,
    const std::string& expected_model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& expected_options
);

} // namespace nova::io
