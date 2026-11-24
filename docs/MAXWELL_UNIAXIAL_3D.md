# Driven 3D uniaxial Maxwell FEM (experimental)

The C++20 API in `nova/physics/maxwell_uniaxial_3d.hpp` extends the driven
edge-element Maxwell problem to cellwise constant reciprocal uniaxial materials.

With the exp(-i omega t) convention, the relative permittivity tensor is

    epsilon = epsilon_o (I - a a^T) + epsilon_e a a^T

where a is the normalized real optic axis. Both principal imaginary parts
must be nonnegative. Principal real parts may be negative. Scalar relative
permeability must be finite and strictly positive. The director is invariant
under nonzero scaling and reversal. The solver expects SI geometry and
vacuum wavenumber in inverse metres.

The discrete operator is curl-curl divided by relative permeability minus
k0^2 times the tensor-weighted edge mass matrix. Four-point degree-two
tetrahedral quadrature exactly integrates the mass of affine first-order
Nedelec basis functions. The assembly shares the existing edge ordering,
PEC elimination and real-block GMRES solve.

The cell absorption diagnostic integrates

    q = (omega epsilon0/2) *
        [Im(epsilon_o) |E_perp|^2 + Im(epsilon_e) |E_parallel|^2]

in watts over each tetrahedron. Both projections use the normalized
director; nonnegative principal losses imply nonnegative absorption.

## Verification cases

The existing Maxwell unit-test translation unit now includes:
- Complex isotropic limit against the existing scalar assembler.
- Reversed/scaled director and permuted element orientation invariance.
- Rotated lossy tensor with an independently prescribed constant field and
  exact manufactured source; exact edge circulations and integrated loss.
- PEC elimination, zero solution and invalid material rejection.

Run the standard CMake build and CTest suite. These tests are not a report
of a completed CI run.

## Limitations

The driven and PEC cavity APIs are limited to uniaxial media; neither is a waveguide eigenmode or port solver. It does not model biaxial/nonreciprocal tensors,
spatially varying directors within a cell, liquid-crystal reorientation,
Kerr nonlinearity, PML or material parameter fitting. Real device predictions
need separate validation and mesh-convergence studies.

## Lossless uniaxial cavity modes

The additional API in `nova/physics/maxwell_cavity_3d.hpp` assembles
the PEC generalized eigenproblem K u = lambda M u, where K contains the
scalar inverse-permeability curl form and M is the uniaxial tensor mass.
The same exact mass quadrature is used for the driven problem.

Both principal permittivities must be real, strictly positive, and
frequency-independent for the physical eigenproblem. Lossy or indefinite
materials must use the driven solver instead. The existing bounded dense
reference eigensolver filters the numerical gradient nullspace, reports
relative residuals, and maps each positive eigenvalue to
f = c sqrt(lambda)/(2 pi).

~~~cpp
#include "nova/physics/maxwell_cavity_3d.hpp"

const auto mesh = nova::fem::make_box_mesh_3d(0,1,0,1,0,1,2,2,2).mesh;
const nova::physics::MaxwellUniaxialMaterial3D medium{
    {2.0,0.0}, {6.0,0.0}, {1.0,2.0,3.0}, 1.0
};
const std::vector<nova::physics::MaxwellUniaxialMaterial3D>
    materials(mesh.elements.size(),medium);
const auto matrices = nova::physics::assemble_maxwell_uniaxial_cavity_3d(
    mesh,materials);
const auto modes = nova::physics::solve_maxwell_cavity_3d(matrices,3);
// Check modes.converged, modes.modes.size(), and residuals before use.
~~~

New automated regression cases check the exact one-free-edge spectrum
against independent barycentric integration, the isotropic limit, the
min-max eigenvalue bounds for principal permittivities 2 and 6, gradient
nullity, the driven K-k0^2 M identity, and rejection of lossy/indefinite
materials. The one-edge benchmark verifies the discrete equations, not
continuum resonator accuracy. Mesh convergence and solver scalability
remain separate validation tasks. No complex, dispersive, open-resonator
or anisotropic waveguide eigenmode support is claimed.
