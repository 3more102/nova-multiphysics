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

} // namespace nova::physics
