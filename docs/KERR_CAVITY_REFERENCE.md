# Single-mode Kerr cavity reference

The header nova/physics/kerr_cavity.hpp implements a dimensionless driven dissipative single-mode cavity with normalized complex amplitude a, time tau in amplitude-decay lifetimes, detuning d in amplitude linewidths, positive normalized Kerr shift, and nonnegative normalized input drive p.

The dynamics are da/dtau = -(1+i(d-|a|^2))a + sqrt(p). The steady intensity x=|a|^2 obeys p=x[1+(d-x)^2]. The phase is recoverable as a=sqrt(p)/(1+i(d-x)). The intensity must lie in [0,p]. The extrema are x=(2d +/- sqrt(d^2-3))/3 when d>sqrt(3). Splitting at these extrema isolates monotonic intervals, followed by bracketed bisection. No initial guess or branch selection is required.

Linearizing the two-real-dimensional dynamics yields Jacobian trace -2 and determinant 1+(d-x)(d-3x). Thus positive differential slope gives linear asymptotic stability, negative gives instability, and zero is nonhyperbolic. The API reports a near-zero slope as marginal using a floating-point tolerance. Near a fold or cusp, roots closer than numerical resolution may merge or their classification may be uncertain; this is not a certified interval-arithmetic solver.

Supported numerical range: |d| <= 1000 and p=0 or 1e-12 <= p <= 1e6. Inputs are normalized, not SI intensities or powers. Nonfinite or out-of-range values are rejected.

Verification: exact cubic factorization for d=3,p=3; exact double root at d=2,p=2; below/above-fold checks; zero-drive and resonant limits; 64 sweep combinations checking root residual, bounds and ordering; invalid input rejection. Run ctest -R nova_kerr_cavity_tests --output-on-failure from the configured build directory.

Scope: a reduced-order reference for future nonlinear frequency-domain optical FEM and photonic resonator continuation, not spatial Maxwell FEM, transient propagation, experimental validation, or multimode stability.
