#pragma once

#include "nova/model/thermal_project.hpp"

#include <filesystem>
#include <iosfwd>

namespace nova::io {

inline constexpr std::size_t max_thermal_project_bytes = 65536;

[[nodiscard]] model::ThermalProject read_thermal_project(std::istream& input);
void write_thermal_project(std::ostream& output, const model::ThermalProject& project);
[[nodiscard]] model::ThermalProject load_thermal_project(const std::filesystem::path& path);
// Publish a validated canonical record atomically on the destination filesystem.
void save_thermal_project(const std::filesystem::path& path, const model::ThermalProject& project);

// Export nodal kelvin, cell Fourier heat flux and volumetric generation.
void save_thermal_result_vtk(const std::filesystem::path& path, const model::ThermalStudyResult& result);

} // namespace nova::io
