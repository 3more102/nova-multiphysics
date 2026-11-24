#pragma once

#include "nova/fem/nedelec_tetrahedron_3d.hpp"
#include "nova/linalg/gmres.hpp"

#include <complex>
#include <span>

namespace nova::physics {

struct MaxwellCellMaterial3D {
    std::complex<double> relative_permittivity{1.0,0.0};
    double relative_permeability{1.0};
};

struct MaxwellSystem3D {
    std::vector<fem::Edge3D> edges; // lexicographic global edge ordering
    std::vector<std::size_t> free_edges; // indices in edges
    linalg::CsrMatrix real_block; // [Re A, -Im A; Im A, Re A]
};

// Time convention exp(-i omega t), k0=omega/c in inverse metres.
// PEC optionally applied on the entire exterior. Otherwise homogeneous
// natural curl boundary conditions are assumed. No ports or PML.
[[nodiscard]] MaxwellSystem3D assemble_maxwell_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials,
    double vacuum_wavenumber, bool pec_on_exterior);

// Homogeneous PEC on selected exterior triangular faces; the remaining
// exterior faces have homogeneous natural curl conditions. Face orientation
// and order are immaterial. Duplicate, interior and nonexistent faces reject.
[[nodiscard]] MaxwellSystem3D assemble_maxwell_3d_selective_pec(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials,
    double vacuum_wavenumber,
    std::span<const fem::BoundaryFace3D> pec_faces);

// Assemble a piecewise constant vector weak source F in
// curl(mu_r^-1 curl E)-k0^2 epsilon_r E=F, in global edge order.
using ComplexVector3D = std::array<std::complex<double>,3>;
[[nodiscard]] std::vector<std::complex<double>> assemble_maxwell_cell_load_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const ComplexVector3D> cell_sources,
    const MaxwellSystem3D& system);

struct MaxwellSolution3D {
    std::vector<std::complex<double>> edge_field; // zero on PEC edges
    linalg::SolverResult diagnostics;
};

// Full-edge weak load is reduced to free DOFs and solved with existing GMRES.
// A nonconverged iterate is returned with diagnostics.converged == false.
[[nodiscard]] MaxwellSolution3D solve_maxwell_3d(
    const MaxwellSystem3D& system,
    std::span<const std::complex<double>> edge_load,
    const linalg::GmresOptions& options = {});


struct MaxwellFieldSample3D {
    ComplexVector3D electric_field{};
    ComplexVector3D magnetic_field{};
    fem::EdgeVector3D time_averaged_poynting{};
    double absorbed_power_density{}; // W/m^3
};

// Reconstruct physical fields on an affine tetrahedron at (r,s,t).
// Coefficients are in the element's local edge order, with global orientation.
[[nodiscard]] MaxwellFieldSample3D evaluate_maxwell_field_3d(
    const fem::NedelecTetrahedron3D& cell,
    const std::array<std::complex<double>,6>& coefficients,
    const MaxwellCellMaterial3D& material, double vacuum_wavenumber,
    double r, double s, double t);

// Integrate dielectric loss exactly for the affine edge field. Returns
// one absorbed power (W) per mesh cell, suitable for conservative transfer.
[[nodiscard]] std::vector<double> integrate_maxwell_absorption_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials,
    double vacuum_wavenumber,
    const MaxwellSystem3D& system,
    std::span<const std::complex<double>> edge_field);

} // namespace nova::physics
