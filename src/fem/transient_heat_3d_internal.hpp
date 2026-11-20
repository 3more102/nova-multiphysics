#pragma once

#include "nova/fem/transient_heat_3d.hpp"
#include "nova/linalg/conjugate_gradient.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace nova::fem::detail {

struct PreparedTransientHeat3D {
    const TetrahedralMesh3D* mesh{};
    ScalarField3DTime time_dependent_source;
    ScalarField3DTime time_dependent_diffusivity;
    SymmetricTensorField3DTime time_dependent_diffusivity_tensor;
    ScalarField3DTime time_dependent_volumetric_capacity;
    std::vector<double> element_jacobians;
    std::vector<bool> constrained;
    std::vector<std::vector<ScalarField3DTime>> time_dependent_prescribed;
    struct TimeDependentNeumannFace {
        BoundaryFace3D face{};
        ScalarField3DTime outward_flux;
        double area_scale{};
    };
    std::vector<TimeDependentNeumannFace> time_dependent_neumann_faces;
    struct TimeDependentRobinFace {
        BoundaryFace3D face{};
        ScalarField3DTime transfer_coefficient;
        ScalarField3DTime ambient_value;
        double area_scale{};
    };
    std::vector<TimeDependentRobinFace> time_dependent_robin_faces;
    linalg::CsrMatrix mass;
    linalg::CsrMatrix stiffness;
    linalg::CsrMatrix robin;
    std::vector<double> volume_load;
    std::vector<double> neumann_load;
    std::vector<double> ambient_load;
    std::vector<double> total_load;
    std::vector<double> capacity_weights;
    std::vector<std::optional<double>> prescribed;
    std::vector<double> initial_state;
};

struct TransientHeat3DMaterialState {
    linalg::CsrMatrix mass;
    // Conductivity stiffness only; static and time-dependent Robin matrices
    // are combined separately when preparing the time-step operators.
    linalg::CsrMatrix stiffness;
    std::vector<double> capacity_weights;
};

struct TransientHeat3DStepOperators {
    linalg::CsrMatrix lhs;
    linalg::CsrMatrix rhs_operator;
    linalg::CsrMatrix system;
    std::vector<linalg::Triplet> boundary_couplings;
};

// One cache is used for full-step trials and another for half-step trials.
// Each slot keeps the most recently prepared dt-dependent/time-level matrix
// and its exact-matrix IC(0) cache. A changed dt or Robin transfer matrix
// rebuilds operators; CachedCG decides whether a new factor is required.
struct TransientHeat3DStepCache {
    std::optional<TransientHeat3DStepOperators> operators;
    double time_step{};
    TransientHeat3DScheme scheme{TransientHeat3DScheme::BackwardEuler};
    linalg::CachedConjugateGradient cached_cg;
    // Dynamic Robin matrices are part of the time-level operators. Exact
    // equality lets ambient-only changes reuse both operators and IC(0).
    std::optional<linalg::CsrMatrix> previous_time_dependent_robin;
    std::optional<linalg::CsrMatrix> next_time_dependent_robin;
    std::optional<linalg::CsrMatrix> previous_material_mass;
    std::optional<linalg::CsrMatrix> next_material_mass;
    std::optional<linalg::CsrMatrix> previous_material_stiffness;
    std::optional<linalg::CsrMatrix> next_material_stiffness;
    std::size_t operator_builds{};
};

[[nodiscard]] PreparedTransientHeat3D prepare_transient_heat_3d(
    const TransientHeat3DProblem& problem,
    double start_time
);

[[nodiscard]] TransientHeat3DResult solve_prepared_transient_heat_3d(
    const PreparedTransientHeat3D& prepared,
    const TransientHeat3DOptions& options,
    std::vector<double> initial_state,
    TransientHeat3DStepCache& cache
);

} // namespace nova::fem::detail
