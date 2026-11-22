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
    // Optional constant applied traction on individually selected exterior faces.
    std::vector<ElasticFaceTraction3D> face_tractions;
};

struct LinearElasticity3DResult {
    std::vector<ElasticVector3D> displacement; // metres, one per vertex
    std::vector<ElasticTensor6> cell_strain;    // dimensionless
    std::vector<ElasticTensor6> cell_stress;    // Pa
    double strain_energy_joules{};             // 0.5 * integral(sigma:epsilon)
    linalg::SolverResult linear_solver;
};

// Experimental 3D continuous P1 tetrahedral, isotropic, linear elasticity.
// Symmetric Galerkin stiffness, consistent constant body/face loads, and
// symmetric displacement elimination. Every connected component must have
// three noncollinear fully prescribed vertices to remove rigid body modes.
// This first version omits thermal strain, geometric/material nonlinearity,
// contact, partial-component constraints, and recovered equilibrated tractions.
[[nodiscard]] LinearElasticity3DResult solve_linear_elasticity_3d(
    const LinearElasticity3DProblem& problem,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::physics
