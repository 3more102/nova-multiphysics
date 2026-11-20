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
    std::vector<double> element_jacobians;
    std::vector<bool> constrained;
    std::vector<std::vector<ScalarField3DTime>> time_dependent_prescribed;
    struct TimeDependentNeumannFace {
        BoundaryFace3D face{};
        ScalarField3DTime outward_flux;
        double area_scale{};
    };
    std::vector<TimeDependentNeumannFace> time_dependent_neumann_faces;
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

struct TransientHeat3DStepOperators {
    linalg::CsrMatrix lhs;
    linalg::CsrMatrix rhs_operator;
    linalg::CsrMatrix system;
    std::vector<linalg::Triplet> boundary_couplings;
};

// One cache is used for full-step trials and another for half-step trials.
// Each slot keeps the most recently prepared dt-dependent matrix and its
// exact-matrix IC(0) cache. A changed dt rebuilds operators; CachedCG decides
// whether the resulting constrained matrix actually requires a new factor.
struct TransientHeat3DStepCache {
    std::optional<TransientHeat3DStepOperators> operators;
    double time_step{};
    TransientHeat3DScheme scheme{TransientHeat3DScheme::BackwardEuler};
    linalg::CachedConjugateGradient cached_cg;
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
