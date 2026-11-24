# Lossless cavity energy and mode volume

The real edge coefficients returned by the PEC cavity solver reconstruct
an affine electric field and a constant curl on each tetrahedron. With
lambda = k0^2 > 0, define E_c = integral(epsilon_r |E|^2) and
H_c = integral(|curl E|^2 / (mu_r lambda)). The physical time-averaged
electric and magnetic energies are epsilon_0/4 times these integrals
when the field is expressed in V/m. Their global equality follows from
the lossless generalized eigenproblem. The reported imbalance is
|sum E_c - sum H_c| / sum(E_c + H_c). It is a diagnostic, not an
independent certificate that a vector is an eigenmode.

The symmetric four-point degree-two quadrature integrates E_c exactly
(up to floating-point error); the constant curl integral is exact.
Cell electric participation is E_c / sum E_c and combined participation
is (E_c + H_c) / sum(E_c + H_c). Both are invariant under amplitude
rescaling. Aggregate cells by material or subdomain to obtain region
participation. The electric effective mode volume is sum E_c divided by
max(epsilon_r |E|^2). Since the squared norm of an affine field is
convex, the maximum on each closed tetrahedron occurs at a vertex;
each adjacent cell's material trace is checked separately.

The API checks lossless positive scalar materials, finite amplitudes,
mesh/edge correspondence, zero constrained coefficients and positive
frequency eigenvalue. Arbitrary nonzero fields are accepted for
diagnostics; the caller must use the matching mesh, material and
converged eigenmode for physical interpretation. This is not a
dispersive, open, lossy, or quasinormal-mode normalization and does not
compute Q factors or Purcell enhancement.

## Regression cases

- Single unit tetrahedron, first oriented edge basis: integral |E|^2 =
  1/12, integral |curl E|^2 = 4/3, peak |E|^2 = 3, volume = 1/36.
  The artificial lambda=16 provides a direct check of both integrals;
  it is not asserted to be a physical PEC eigenmode.
- Heterogeneous positive dielectric cube: M-normalization, electric/
  magnetic equipartition, participation sums, and amplitude invariance.
- Invalid loss, zero lambda and zero field are rejected.

Run the existing CMake/CTest suite. No additional dependencies.
