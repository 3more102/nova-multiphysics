#pragma once

#include "nova/fem/mesh3d.hpp"
#include "nova/linalg/conjugate_gradient.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace nova::physics {

using ElasticVector3D = std::array<double, 3>;
using ElasticTensor6 = std::array<double, 6>; // xx, yy, zz, xy, yz, xz; tensor shear

struct IsotropicElasticMaterial3D {
    double young_modulus{};  // Pa
    double poisson_ratio{};  // -1 < nu < 0.5
};

struct DisplacementConstraint3D {
    std::size_t node{};
    ElasticVector3D value{}; // prescribe all three Cartesian components, metres
};

struct DisplacementComponent3D {
    std::size_t node{};
    std::size_t component{}; // 0=x, 1=y, 2=z in global Cartesian axes
    double value{};          // metres
};

struct ElasticFaceTraction3D {
    fem::BoundaryFace3D face{};
    ElasticVector3D traction{}; // applied global Cartesian traction, Pa
};

struct LinearElasticity3DProblem {
    fem::TetrahedralMesh3D mesh;
    // Exactly one isotropic, small-strain linear material per tetrahedron.
    std::vector<IsotropicElasticMaterial3D> cell_materials;
    // Fully constrain all three DOFs at each selected vertex.
    std::vector<DisplacementConstraint3D> prescribed_displacements;
    // Optional one constant body-force density (N/m^3) per tetrahedron.
    std::vector<ElasticVector3D> cell_body_force_density;
    // Optional isotropic thermal eigenstrain: alpha[e] * delta_T[e].
    // Supply both cellwise arrays (one entry per tetrahedron), or neither.
    // Negative temperature changes and negative expansion coefficients are valid.
    std::vector<double> cell_thermal_expansion_per_kelvin;
    std::vector<double> cell_temperature_change_kelvin;
    // Optional constant applied traction on individually selected exterior faces.
    std::vector<ElasticFaceTraction3D> face_tractions;
    std::vector<DisplacementComponent3D> component_constraints{};
};

struct LinearElasticity3DResult {
    std::vector<ElasticVector3D> displacement; // metres, one per vertex
    std::vector<ElasticTensor6> cell_strain;    // total symmetric displacement strain
    std::vector<ElasticTensor6> cell_stress;    // C:(cell_strain - thermal_strain), Pa
    double strain_energy_joules{};             // 0.5 * integral(sigma:elastic_strain)
    linalg::SolverResult linear_solver;
    std::vector<ElasticVector3D> support_reaction; // internal minus applied force [N], zero on free components
};

// Principal stresses in descending algebraic order and von Mises equivalent stress.
[[nodiscard]] ElasticVector3D principal_stresses_3d(const ElasticTensor6& stress);
[[nodiscard]] double von_mises_stress_3d(const ElasticTensor6& stress);

// Experimental 3D continuous P1 tetrahedral, isotropic, linear elasticity.
// Symmetric Galerkin stiffness, consistent constant body/face loads, and
// symmetric displacement elimination. Supports must suppress all six rigid
// modes on each face-connected solid.
// Includes prescribed, cellwise isotropic thermal eigenstrain; temperature
// is an input, not a solved or bidirectionally coupled thermal field.
// Omits geometric/material nonlinearity, contact, arbitrary-direction constraints, and recovered equilibrated tractions.
[[nodiscard]] LinearElasticity3DResult solve_linear_elasticity_3d(
    const LinearElasticity3DProblem& problem,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::physics
