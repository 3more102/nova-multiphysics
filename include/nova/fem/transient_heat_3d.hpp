#pragma once

#include "nova/fem/conductivity_3d.hpp"
#include "nova/fem/h1_space_3d.hpp"
#include "nova/fem/transient_heat.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

using TransientScalarField3D =
    std::function<double(double, double, double, double)>;

// Consistent Galerkin mass and stiffness matrices for continuous P1 tetrahedra.
using ScalarH1Operators3D = ScalarH1Operators;

[[nodiscard]] ScalarH1Operators3D assemble_scalar_h1_operators_3d(
    const TetrahedralMesh3D& mesh,
    const SpatialScalarField3D& volumetric_capacity,
    const SpatialScalarField3D& conductivity
);

// Assemble the consistent mass and stiffness for a spatial SPD tensor field.
[[nodiscard]] ScalarH1Operators3D assemble_anisotropic_h1_operators_3d(
    const TetrahedralMesh3D& mesh,
    const SpatialScalarField3D& volumetric_capacity,
    const SymmetricTensorField3D& conductivity_tensor
);

struct TransientDirichletFace3D {
    BoundaryFace3D face{};
    TransientScalarField3D temperature;
};

struct TransientNeumannFace3D {
    BoundaryFace3D face{};
    // Outward conductive heat flux q_n = -(K grad(T)) dot n.
    TransientScalarField3D outward_flux;
};

struct TransientRobinFace3D {
    BoundaryFace3D face{};
    // Outward flux = h * (T - ambient), with finite h >= 0.
    TransientScalarField3D transfer_coefficient;
    TransientScalarField3D ambient_temperature;
};

struct TransientHeat3DProblem {
    TetrahedralMesh3D mesh;
    // Volumetric heat capacity rho*c > 0 and conductivity k > 0 are spatial,
    // time-independent fields.
    SpatialScalarField3D volumetric_capacity;
    SpatialScalarField3D conductivity;
    TransientScalarField3D source;
    SpatialScalarField3D initial_temperature;
    // Optional exact P1 nodal restart state, overriding initial_temperature.
    std::vector<double> initial_dof_values;
    // Boundary faces must be exterior and mutually exclusive. Other faces are
    // insulated (zero outward flux).
    std::vector<TransientDirichletFace3D> dirichlet_faces;
    std::vector<TransientNeumannFace3D> neumann_faces;
    std::vector<TransientRobinFace3D> robin_faces;
    // Supply exactly one of conductivity (k I) and conductivity_tensor (K).
    // Both are time independent and sampled only when assembling the operators.
    SymmetricTensorField3D conductivity_tensor{};
};

using TransientHeat3DOptions = TransientHeatTimeOptions;
using TransientHeatStepDiagnostics3D = TransientHeatStepDiagnostics;
using TransientHeat3DResult = TransientHeatResult;

// Solves rho*c*dT/dt - div(K grad(T)) = Q using tetrahedral P1 elements,
// Backward Euler or Crank-Nicolson, and time-dependent exterior face data.
[[nodiscard]] TransientHeat3DResult solve_transient_heat_3d(
    const TransientHeat3DProblem& problem,
    const TransientHeat3DOptions& options = {}
);

} // namespace nova::fem
