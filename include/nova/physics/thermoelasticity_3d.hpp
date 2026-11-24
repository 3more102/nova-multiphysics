#pragma once

#include "nova/physics/linear_elasticity_3d.hpp"

#include <memory>
#include <span>

namespace nova::physics {

struct ThermalExpansion3D {
    double per_kelvin{};
    double reference_temperature{}; // stress-free temperature, same scale as the nodal field
};

// Monotonic wall-clock seconds. Setup is measured once. The final three
// counters accumulate completed response calls, including returned CG failures.
struct ElasticityTimings3D {
    double topology{};
    double element_setup{};       // material, gradient, constraint and static load preparation
    double stiffness_assembly{};  // element blocks and global constrained CSR, fused by row
    double thermal_load{};        // field transfer, RHS and warm-start tolerance
    double linear_solve{};        // CG including preconditioner work
    double recovery{};            // displacement, stress/strain, energy and warm-state copy
};

// Quasistatic P1 mechanics on a fixed mesh. The object owns its assembled
// geometry, loads, constraints and CSR stiffness. Each solve updates only the
// eigenstrain load and recovery, warm-starting from the last converged field.
// Relative tolerance for a repeated thermal solve uses the cold-start free-force
// norm, so an unchanged load has the same accuracy target. Use one instance per
// concurrent study.
class PreparedElasticity3D {
public:
    PreparedElasticity3D(
        const fem::TetrahedralMesh3D& mesh,
        std::span<const IsotropicElasticMaterial3D> cell_materials,
        std::span<const DisplacementConstraint3D> constraints,
        std::span<const ElasticVector3D> cell_body_forces = {},
        std::span<const ElasticFaceTraction3D> tractions = {},
        std::span<const DisplacementComponent3D> component_constraints = {}
    );
    ~PreparedElasticity3D();
    PreparedElasticity3D(PreparedElasticity3D&&) noexcept;
    PreparedElasticity3D& operator=(PreparedElasticity3D&&) noexcept;

    // Temperature has one DOF per mesh vertex. Stress is the exact cell mean
    // of the discrete stress. Energy includes within-cell thermal variation.
    [[nodiscard]] LinearElasticity3DResult solve_temperature(
        std::span<const double> temperature,
        std::span<const ThermalExpansion3D> cell_expansion,
        const linalg::SolverOptions& options = {}
    );

    // A prescribed, constant isotropic eigenstrain per cell. Empty means zero.
    [[nodiscard]] LinearElasticity3DResult solve_cell_eigenstrain(
        std::span<const double> eigenstrain = {},
        const linalg::SolverOptions& options = {}
    );

    [[nodiscard]] ElasticityTimings3D timings() const noexcept;
    [[nodiscard]] std::size_t stiffness_nonzeros() const noexcept;
    [[nodiscard]] std::size_t factor_builds() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace nova::physics
