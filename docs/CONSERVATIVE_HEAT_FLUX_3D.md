# Conservative steady 3D face-integrated heat-flow reconstruction

Status: experimental independent postprocessor, not a replacement for the P1
conduction solve, the weak-form global audit, or a mixed FEM method.

## Purpose and definitions

The raw P1 tetrahedron flux, q_K = -K grad(T_h), generally does not have
continuous normal components across interior faces and generally does not
satisfy elementwise source balance. This postprocessor recovers one integrated
heat rate Q_f in W per unique face, with positive sign outward from its stored
owner cell. On an interior face the neighboring element sees exactly -Q_f.

The reconstruction enforces for every tetrahedron K:

    sum_{f on boundary(K)} Q_{K,f} = integral_K s(x) dx.

The source integral uses the same positive degree-two tetrahedral quadrature
as the stationary FEM assembly, including per-cell source additions.

## What is fixed, and what can change?

* Each Neumann exterior face is fixed at the quadrature-integrated prescribed
  outward heat rate.
* Each Robin exterior face is fixed at integral h(x)(T_h(x)-T_ambient(x)) dA,
  evaluated using the supplied P1 nodal field.
* Unassigned exterior faces are fixed at zero (insulated).
* Explicitly tagged Dirichlet exterior faces are free heat outlets.
* Each interior face is a single adjustable heat-transfer rate, shared by
  the two adjacent elements with opposite signs.

A vertex-only Dirichlet constraint is **not** an exterior-face label. If energy
is removed at an interior pinned vertex, a face-only conservation model cannot
represent that sink without extra source/reaction geometry. The API therefore
rejects an unbalanced closed component, rather than inventing face rates.

## Reconstruction method

For a face with initial outward rate Q0 (towards the neighbor if interior),
the target is the mean of the two cells' raw outward-normal rates, expressed
in the owner's orientation. Dirichlet-face targets use the owner's raw normal
rate. Fixed boundary rates are used as-is.

Find the free rate vector that minimizes:

    0.5 * sum_{f free} (Q_f - Q0_f)^2,

subject to all elementwise balances. With cell potentials p_K, the correction
has the form

    Q_f = Q0_f + p_owner - p_neighbor   (interior)
    Q_f = Q0_f + p_owner                (Dirichlet exterior).

This yields a graph-Laplacian linear system. We solve it with existing CSR
Conjugate Gradient and Jacobi preconditioning. A component with no Dirichlet
face has a constant-potential nullspace and requires net prescribed heat
balance; one potential is then gauge-pinned without introducing heat exchange.

The choice of equal weights makes the constrained least-squares problem
well-defined, but does not constitute an accuracy-optimal recovery estimate.
There is no pointwise H(div) vector field, normal-density interpolation,
guaranteed error bound, or time-dependent capacity term in this API.

## Public interface

Include `nova/fem/conservative_heat_flux_3d.hpp` and call
`reconstruct_conservative_heat_flux_3d(problem, temperature)`.

The result provides:

* unique faces with canonical vertex IDs, owner/neighbor cell indices, and
  owner-outward integrated heat rates in W;
* integrated generated and reconstructed outward heat per cell;
* maximum absolute element balance defect and correction iteration count.

A shared face appears exactly once; add its rate to the owner balance and
subtract it from the neighbor balance. All cases must pass both the existing
FE weak-form validation and this reconstruction's compatibility checks.

## Verification scope

The test suite covers single-cell volumetric generation with a tagged
Dirichlet outlet, pure Neumann compatibility and prescribed flux preservation,
incompatible closed volumes, orientation invariance, anisotropic affine
patch behavior, and a multi-cell anisotropic two-material heated wall with
Robin convection. These are regression checks, not a substitute for
independent cross-code benchmark verification.

Follow-on work: area- or material-weighted objective, precise physical
Dirichlet reaction localization, H(div)-conforming RT0/BDM recovery, flux-jump
a posteriori estimators, and transient storage-rate reconstruction.
