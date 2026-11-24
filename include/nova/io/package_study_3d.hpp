#pragma once
#include "nova/physics/package_3d.hpp"
#include <filesystem>
#include <iosfwd>
#include <string>

namespace nova::io {
struct PackageProbe3D {
    std::string name;
    std::size_t cell{};
    std::array<double,4> weights{};
};
struct PackageStudy3D {
    physics::PackageModel3D model;
    physics::PackageOptions3D options;
    std::vector<int> material_tags, terminal_tags; // cells, then prescribed voltage nodes
    std::vector<PackageProbe3D> probes;
    double initial_kelvin{}, end_time{};
    std::size_t frame_stride{1};
    std::string model_id;
};
[[nodiscard]] PackageStudy3D read_package_study_3d(const std::filesystem::path& model_file);

// The exact mesh and model identity, two accepted temperature levels, algebraic
// potential and next proposed step are verified before a state can be restored.
void write_package_checkpoint_3d(std::ostream& out, const fem::TetrahedralMesh3D& mesh,
    const std::string& model_id, const physics::PackageState3D& state);
[[nodiscard]] physics::PackageState3D read_package_checkpoint_3d(std::istream& in,
    const fem::TetrahedralMesh3D& mesh,const std::string& model_id);
}
