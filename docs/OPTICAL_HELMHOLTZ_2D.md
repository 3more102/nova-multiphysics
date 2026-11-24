# Experimental scalar complex Helmholtz FEM in 2D

This small reference path is independent of the planar 1D optical solver. It uses
the existing validated conforming triangular mesh, P1 nodal fields, exact
constant-coefficient stiffness and consistent mass integrals, and a bounded
dense complex solver with row equilibration and partial pivoting. Maximum
1024 nodes; O(N^2) storage and O(N^3) worst-case work. Not a production solver.

For z-invariant isotropic scalar materials and exp(-i omega t), the out-of-plane
electric (TE_Ez) and magnetic (TM_Hz) reductions solve
-div(a grad u) - k0^2 b u = f, where (a,b)=(1/mu_r,epsilon_r) for Ez and
(1/epsilon_r,mu_r) for Hz. The names specify the actual scalar, avoiding
ambiguity with planar plane-of-incidence naming. k0 is positive in rad/m,
coordinates in metres. Passive epsilon and mu have nonnegative imaginary
parts. Exact zero inverse-weight coefficients are rejected; zero mass-weight
coefficients are permitted. Conductivity may be represented by positive
imaginary epsilon_r=sigma/(omega epsilon0) plus dielectric contribution.

For conjugated test function v, the weak equation is
integral(a grad u dot grad conj(v) - k0^2 b u conj(v))
+ integral_boundary(alpha u conj(v))
= integral(f conj(v)) + integral_boundary(g conj(v)).
Boundary convention a du/dn + alpha u = g; unassigned edges have zero weighted
normal derivative. Exterior nodal Dirichlet data are strongly eliminated.
For normally outgoing homogeneous waves alpha=-i*a*q where
q=k0*sqrt(epsilon_r*mu_r) is the appropriate outgoing branch; at a left
incident plane g=-2*i*a*q for unit incident scalar amplitude.
This local impedance is not a general exact oblique radiation condition.
Piecewise constant coefficients on a conforming interface enforce field and
weighted normal derivative continuity in the weak sense.

The result contains nodal scalar field, constant cell gradients, free-row
absolute infinity residual and normwise backward error
||Au-b||inf/(||A||inf ||u||inf+||b||inf) on the original free rows including
prescribed columns. A small backward error is not a forward error guarantee.

Cell loss is the exact P1 integral of -Im(a)|grad u|^2 +
k0^2 Im(b)|u|^2, not watts. For physical Ez, division by 2 omega mu0 gives
W/m; for physical Hz use 2 omega epsilon0. No incident power normalization
or PML loss separation is claimed. No anisotropy, full-vector coupling,
general ports, eigenmodes, PML, or large-mesh solver is implemented here.

Regression cases: normally propagating homogeneous waves with analytical
exp(i*q*x), both scalar formulations, two mesh sizes, reversed orientation,
electric and magnetic dissipation against the exponential integral, exact
affine Ez at epsilon=0 and singular Hz rejection, invalid interior edges and
active material rejection. Run cmake --build <build> then ctest --test-dir
<build> --output-on-failure. Full CI and numerical results pending.
