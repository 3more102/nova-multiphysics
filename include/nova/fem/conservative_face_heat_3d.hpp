#pragma once

#include "nova/fem/diffusion_3d.hpp"
#include "nova/linalg/solver_types.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace nova::fem {

// One signed integrated heat rate per geometrical face, in watts.
// For an interior face it is positive leaving first_cell and entering
// second_cell. No second_cell is present for an exterior face.
enum class FaceHeatCondition3D {
    Interior, Dirichlet, Neumann, Robin, Insulated
};

struct ReconciledFaceHeat3D {
    BoundaryFace3D nodes{}; // sorted vertex indices
    std::size_t first_cell{};
    std::size_t second_cell{}; // meaningful only if interior == true
    bool interior{};
    FaceHeatCondition3D condition{FaceHeatCondition3D::Insulated};
    double area{};
    double outward_heat{};     // final single-valued heat rate [W]
    double raw_outward_heat{}; // cell-averaged P1 estimate before correction [W]
};

struct FaceHeatReconciliation3D {
    std::vector<ReconciledFaceHeat3D> faces;
    std::vector<double> generated_heat_per_cell; // [W]
    std::vector<double> cell_balance_defect;     // generation - outward [W]
    double max_absolute_cell_defect{};           // [W]
    double max_absolute_correction{};            // [W]
    double correction_l2_flux_density{};         // sqrt(sum_F area_F*(delta_F/area_F)^2)
    double total_generated_heat{};               // [W]
    double total_exterior_outward_heat{};         // [W]
    std::size_t floating_components{};
    linalg::SolverResult multiplier_solver;
};

// Constructs conservative FACE-INTEGRATED rates by a minimum-correction
// constrained least-squares problem on the tetrahedral dual graph.
//
// Unknowns are interior face rates and Dirichlet-boundary face rates.
// Explicit Neumann and Robin integrated rates remain fixed; unassigned
// exterior faces remain insulated (zero). Correction minimizes
// sum_F ((new_rate - raw_rate) / face_area)^2, subject to each tetrahedron's
// integrated heat generation equaling its signed outgoing face rates.
//
// The raw interior rate is the average of the two opposite P1 cell rates.
// Output is a single integrated normal rate per face, NOT an H(div)-conforming
// flux field, a pointwise conservation proof, or a flux error bound. The
// least-squares rate is not unique without the stated minimum-norm choice.
//
// Pure-Neumann connected components must have physically balanced source
// and boundary rates. Explicit *nodal* Dirichlet constraints are rejected:
// their concentrated nodal reactions cannot be assigned to exterior faces.
// Use exterior Dirichlet face conditions for this diagnostic instead.
//
// Reuses NOVA's validated mesh, material, volume/face quadrature and P1 flux
// recovery. Does not change or invoke the diffusion solver.
[[nodiscard]] FaceHeatReconciliation3D reconcile_p1_face_heat_rates_3d(
    const ScalarDiffusion3DProblem& problem,
    std::span<const double> nodal_values,
    const linalg::SolverOptions& correction_solver_options = {}
);

} // namespace nova::fem
