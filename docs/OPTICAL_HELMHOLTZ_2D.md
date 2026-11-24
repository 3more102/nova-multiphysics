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

Independent SciPy reference (separate physical-coordinate polynomial
interpolation, three-point triangle quadrature, dense LAPACK solve):
lossless Ez/Hz maximum nodal errors 0.01105535 (52 nodes) and 0.00314519
(175 nodes), refinement ratio 0.284495; lossy Ez/Hz error 0.00151875
(297 nodes). Normalized lossy Ez integral 0.323999200 versus 0.324262469
analytical, Hz 0.144149835 versus 0.144266623 analytical; relative errors
0.0812% and 0.0810%. Independent backward errors below 2.5e-16.
An initial x-only refinement failed the intended ratio (0.4398), so the
regression now refines both axes rather than weakening the tolerance.
These are independent reference calculations, not execution of NOVA C++.

Additional independently evaluated cases: oblique wave exp(i(1.2*x+1.6*y))
on a deterministically perturbed 16x16 triangular grid (289 nodes) had
maximum nodal error 0.00145286. A dielectric step at x=0.5 with epsilon
1 and 4, mu=1, on a 40x10 grid (451 nodes), gave Ez error 0.00138967 and
Hz error 0.00280151 against the exact piecewise plane waves. With Y=a*q,
scalar interface r=(Y1-Y2)/(Y1+Y2), t=1+r. The analytical power check is
r^2+(Y2/Y1)t^2=1 for real positive admittances. This check does not yet
extract numerical S-parameters. A constant complex field with volume source
f=-k0^2*b*u is also registered as an exact manufactured test. C++ execution
of these additions remains pending.
