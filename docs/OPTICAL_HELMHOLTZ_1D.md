# One-dimensional optical Helmholtz FEM

The C++20 API in nova/physics/optical_helmholtz_1d.hpp solves normal-incidence scattering on a piecewise-constant passive, isotropic, nonmagnetic material mesh. It returns complex nodal E fields, complex reflection and transmission amplitudes, power fractions, individual element absorptance and a signed global energy residual.

With s=k0*x and exp(-i omega t), the dimensionless equation is -E''-epsilon_r E=0. The consistent P1 element matrix has diagonal 1/h-epsilon_r*h/3 and off-diagonal -1/h-epsilon_r*h/6. Outgoing impedance conditions are d_n E-i*n*E=0 at the exit and d_n E-i*n*E=-2i*n at the entrance for unit incoming E. The endpoint terms subtract i*n from the diagonal. A complex tridiagonal adjacent-pivoted direct solver uses O(N) storage and arithmetic. Reference planes are the mesh endpoints; exterior indices must be positive and real.

R=|E_left-1|^2, T=(n_exit/n_incident)|E_right|^2, A=sum Im(epsilon_r)*integral |E_h|^2 ds/n_incident. Element integrals are exact for the P1 field. The reported residual 1-R-T-A is computed independently of the absorption integral, not used to define A. A small residual verifies discrete power balance but does not imply spatial accuracy; refine the mesh and compare fields and powers. The solver does not guarantee accuracy on under-resolved or ill-conditioned meshes.

CTest target nova_optical_helmholtz_1d_tests compares complex amplitudes and absorptance against the existing analytical layer cascade for transparent, absorbing and metal-like slabs, checks uniform wave convergence, power conservation, pivoting and invalid input. Run with the existing default build and ctest presets.

Scope: driven 1D scalar FEM, not a general 2D/3D vector Maxwell solver, eigenmode solver, PML, or self-consistent nonlinear solver.

## Graded multilayer verification

The regression suite additionally exercises three dielectric/metal-like layers with unequal thicknesses, unequal real exterior indices (1.3 and 1.8), and a nonuniform mesh with local node coordinate proportional to (j/m)^1.25. An independent exact solution propagates the state [E,dE/ds] backward across each homogeneous layer using complex sine and cosine; the exit outgoing wave is normalized to unit incoming electric amplitude at the entrance. This moderate-thickness test reference is not a replacement for a stable scattering cascade in thick, highly attenuating stacks.

Checks cover complex scattering amplitudes, second-order transmission convergence between 80 and 160 elements per layer, electric fields at material interfaces, individual layer absorption against the exact difference of normalized Poynting flux, global discrete energy balance, translation of the reference planes with the mesh, and reciprocal power transmission under reversal. The reversed mesh is regraded rather than geometrically mirrored, so its reciprocity tolerance includes discretization error. All checks run in the existing Helmholtz CTest target.

## Repeated-pivot algebra verification

A separate regression uses deliberately coarse 13-node lossless and lossy meshes to force repeated row exchanges. It independently assembles a dense weak-form matrix using two-point Gauss quadrature and solves with full-row Gaussian elimination and partial pivoting. It checks that at least two row exchanges occur, compares every nodal field and both complex scattering amplitudes, evaluates the original assembled equation residual, and checks discrete power balance. This is an algebraic stress test, not a claim of continuum accuracy on coarse meshes.
