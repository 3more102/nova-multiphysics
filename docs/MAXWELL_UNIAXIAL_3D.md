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

This is a driven frequency-domain material extension, not an anisotropic
eigenmode or port solver. It does not model biaxial/nonreciprocal tensors,
spatially varying directors within a cell, liquid-crystal reorientation,
Kerr nonlinearity, PML or material parameter fitting. Real device predictions
need separate validation and mesh-convergence studies.
