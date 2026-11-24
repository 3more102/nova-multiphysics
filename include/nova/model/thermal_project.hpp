#pragma once

#include "nova/fem/conductivity_3d.hpp"
#include "nova/fem/energy_balance_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"

#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace nova::model {

// These geometry selections survive changes to box dimensions and mesh density.
enum class BoxSide { XMin, XMax, YMin, YMax, ZMin, ZMax };

struct BoxGeometry {
    std::string id{"box"};
    std::array<double, 3> lengths_m{1.0, 1.0, 1.0};
};

struct StructuredTetMesh {
    std::string id{"mesh"};
    std::array<std::size_t, 3> divisions{1, 1, 1};
};

struct ThermalMaterial {
    std::string id{"solid"};
    // Used when conductivity_tensor_w_m_k is empty.
    double conductivity_w_m_k{1.0};
    // Optional homogeneous symmetric tensor in global Cartesian axes.
    // When present, the scalar coefficient is not used by the solver.
    std::optional<fem::SymmetricTensor3D> conductivity_tensor_w_m_k{};
};

// Preserve the original SDK type name for source compatibility.
using IsotropicThermalMaterial = ThermalMaterial;

struct FixedTemperature { double kelvin; };
struct OutwardHeatFlux { double watts_per_m2; };
struct Convection { double film_w_m2_k; double ambient_kelvin; };
using ThermalBoundary = std::variant<FixedTemperature, OutwardHeatFlux, Convection>;

// Schemas 1 and 2 contain one origin-aligned box, mesh, homogeneous material and stationary study.
// Missing side assignments mean zero outward heat flux.
struct ThermalProject {
    std::string id{"project"};
    BoxGeometry geometry;
    StructuredTetMesh mesh;
    ThermalMaterial material;
    std::string study_id{"heat"};
    double generation_w_m3{};
    std::map<BoxSide, ThermalBoundary> boundaries;
    linalg::SolverOptions solver{1000, 1e-10, 1e-14, linalg::PreconditionerType::Jacobi};
};

inline constexpr std::size_t max_project_tetrahedra = 196608;
inline constexpr std::size_t max_project_iterations = 100000;

void validate_thermal_project(const ThermalProject& project);
[[nodiscard]] fem::ScalarDiffusion3DProblem prepare_thermal_project(const ThermalProject& project);

struct ThermalStudyResult {
    fem::TetrahedralMesh3D mesh;
    fem::ScalarDiffusion3DResult temperature;
    std::vector<fem::P1TetrahedralFlux3D> cell_fluxes;
    fem::DiffusionEnergyBalance3D balance;
    double generation_w_m3{};
};

// Throws if the model is invalid, CG fails to converge, or a nodal absolute
// temperature is outside its physical range.
[[nodiscard]] ThermalStudyResult solve_thermal_project(const ThermalProject& project);

} // namespace nova::model
