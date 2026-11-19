#pragma once

#include "nova/fem/diffusion_3d.hpp"

#include <vector>

namespace nova::physics {

// Experimental electrostatic potential on conforming P1 tetrahedral meshes:
//   -div(epsilon grad(V)) = rho_free, E = -grad(V), D = epsilon E.
// SI units: V [volt], epsilon [F/m], rho_free [C/m^3], D [C/m^2].
// Unassigned exterior faces have zero normal electric displacement.
struct Electrostatics3DProblem {
    fem::TetrahedralMesh3D mesh;

    // Exactly one FEM material representation is required. Each scalar must
    // be positive and each symmetric tensor positive definite.
    fem::ScalarField3D permittivity;
    fem::SymmetricTensorField3D permittivity_tensor;
    std::vector<double> cell_permittivity;
    std::vector<fem::SymmetricTensor3D> cell_permittivity_tensors;

    // Missing spatial charge means zero. Optional cell charge is additive.
    fem::ScalarField3D free_charge_density;
    std::vector<double> cell_free_charge_density;

    std::vector<fem::DirichletConstraint3D> potential_constraints;
    std::vector<fem::DirichletBoundaryFace3D> potential_faces;
    // D dot outward_normal = prescribed surface displacement [C/m^2].
    std::vector<fem::NeumannBoundaryFace3D> displacement_flux_faces;
};

struct Electrostatics3DResult {
    std::vector<double> potential;             // At mesh vertices [V].
    std::vector<fem::Point3D> electric_field;  // Per cell, -grad(V) [V/m].
    std::vector<fem::Point3D> displacement;    // Per cell, averaged D [C/m^2].
    double stored_energy_joules{};             // 1/2 integral E dot D dV.
    linalg::SolverResult linear_solver;
};

// Reuses the validated scalar/tensor P1 Poisson solver and flux quadrature;
// electric displacement is the FEM element-averaged constitutive flux, not
// an equilibrated numerical face flux. No nonlinear dielectric, electrodes
// with floating unknown potentials or AC/time-dependent fields yet.
[[nodiscard]] Electrostatics3DResult solve_electrostatics_3d(
    const Electrostatics3DProblem& problem,
    const linalg::SolverOptions& options = {}
);

} // namespace nova::physics
