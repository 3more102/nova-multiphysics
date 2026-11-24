#pragma once

#include "nova/fem/p1_flux_3d.hpp"
#include "nova/nonlinear/newton_solver.hpp"
#include "nova/physics/thermoelasticity_3d.hpp"
#include "nova/time/implicit_formula.hpp"

#include <functional>
#include <memory>
#include <optional>
#include <span>

namespace nova::physics {

// p(T) = value * exp(slope * (T - reference_temperature)); SI units.
// Electrical value zero selects an electrically insulating cell exactly.
struct ExponentialProperty3D { double value{}, slope{}; };
struct PackageMaterial3D {
    double reference_temperature{300};
    ExponentialProperty3D electrical, thermal, volumetric_capacity, young;
    double poisson{}, expansion{};
    // Optional dimensionless SPD tensors in global axes, multiplied by p(T).
    std::optional<fem::SymmetricTensor3D> electrical_shape{}, thermal_shape{};
};
using TimeSignal3D = std::function<double(double)>;
struct NodalSignal3D { std::size_t node{}; TimeSignal3D value; };
struct PackageSurface3D {
    fem::BoundaryFace3D face{};
    double convection{}, emissivity{};
    TimeSignal3D ambient_kelvin;
    TimeSignal3D outward_flux; // additive prescribed W/m^2
};
struct PackageModel3D {
    fem::TetrahedralMesh3D mesh;
    std::vector<PackageMaterial3D> materials;
    std::vector<std::size_t> cell_material;
    std::vector<NodalSignal3D> voltage, temperature;
    std::vector<PackageSurface3D> cooling;
    // Optional additive uniform source per material, W/m^3, indexed like materials.
    std::vector<TimeSignal3D> heat_sources;
    std::vector<DisplacementConstraint3D> clamps;
    std::vector<DisplacementComponent3D> supports;
    std::vector<ElasticVector3D> body_forces;
    std::vector<ElasticFaceTraction3D> tractions;
};
struct PackageState3D {
    double time{}, previous_time{}, next_dt{};
    std::vector<double> temperature, voltage, previous_temperature;
};
struct TransportFields3D {
    std::vector<fem::P1TetrahedralFlux3D> heat_flux;
    std::vector<fem::Point3D> current_density;
    std::vector<double> joule_density, conductivity, young_modulus;
    std::vector<unsigned char> electrically_active;
    double joule_power{}, source_power{}, cooling_power{}, enthalpy{};
    double terminal_power{}, net_terminal_current{};
    std::vector<double> terminal_current; // same order as model.voltage, positive into domain
};
struct PackageStep3D {
    std::size_t iterations{}, rejected_before{};
    double residual{}, normalized_error{}, energy_defect_watts{}, dirichlet_heat_watts{};
    TransportFields3D transport;
    LinearElasticity3DResult mechanics;
};
enum class TransportStrategy3D { MonolithicNewton, Segregated };
struct PackageOptions3D {
    time::ImplicitScheme scheme{time::ImplicitScheme::BDF2};
    TransportStrategy3D strategy{TransportStrategy3D::MonolithicNewton};
    double initial_dt{0.01}, minimum_dt{1e-7}, maximum_dt{0.1};
    // Zero disables embedded error control. Nonlinear failure can still cut the step.
    double error_kelvin{}, relative_error{};
    double temperature_scale{1}, voltage_scale{1};
    double relative_residual{}, relaxation{1};
    nonlinear::NewtonOptions newton{};
    linalg::SolverOptions mechanics{3000,1e-10,1e-14,linalg::PreconditionerType::Jacobi};
    std::size_t maximum_attempts{10000};
    std::vector<double> events;
};
struct TransportAssembly3D {
    std::vector<double> residual;         // T rows (W), V rows (A), then essential rows replaced
    linalg::CsrMatrix jacobian;           // empty unless requested
    double energy_defect_watts{}, dirichlet_heat_watts{};
};

// One fixed P1 mesh, all fields indexed by original vertex/element order.
// Evaluation uses [T_0...T_n, V_0...V_n]. Inactive potentials have identity rows.
class PreparedPackage3D {
public:
    explicit PreparedPackage3D(PackageModel3D model);
    ~PreparedPackage3D();
    PreparedPackage3D(PreparedPackage3D&&) noexcept;
    PreparedPackage3D& operator=(PreparedPackage3D&&) noexcept;
    [[nodiscard]] const PackageModel3D& model() const noexcept;

    [[nodiscard]] TransportAssembly3D assemble(std::span<const double> unknown,
        const PackageState3D& old, double new_time, time::ImplicitFormula formula, bool tangent) const;
    [[nodiscard]] TransportFields3D recover(std::span<const double> temperature,
        std::span<const double> voltage, double time) const;
    // Consistent algebraic voltage and prescribed temperature at t0.
    [[nodiscard]] PackageState3D initialize(std::vector<double> temperature, double t0 = 0) const;
    // Each complete transport + mechanics equilibrium is committed before calling the observer.
    // If recovery is exhausted, throws with state still at the last accepted time.
    void advance(PackageState3D& state, double end_time, const PackageOptions3D& options,
        const std::function<void(const PackageState3D&, const PackageStep3D&)>& accepted = {});
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
