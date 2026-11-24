# Discrete electromagnetic power audit

For peak phasors with exp(-i omega t), let A e = b denote the reduced
complex Maxwell FEM system, with b assembled from F = i omega mu0 J.
With homogeneous PEC or homogeneous natural curl boundary conditions,
there is no boundary power term. Delivered source power is
-Im(e^H b)/(2 omega mu0). Operator absorbed power is
-Im(e^H A e)/(2 omega mu0). Their difference is the work of the algebraic
residual b-Ae. Cell absorbed power is independently integrated using
the existing positive tetrahedral quadrature.

The API audit_maxwell_power_3d returns all three powers (W), residual
work (W), the Euclidean algebraic residual norm, and two dimensionless
relative discrepancies. The latter use max(abs(a),abs(b)) as denominator
and return zero if both quantities vanish. Near-zero or lossless cases
should be judged with absolute power tolerances rather than relative
discrepancies. A small power residual alone does not establish convergence:
the residual may be orthogonal to the field. Inspect the residual norm
and the underlying GMRES diagnostics as well.

The regression uses a heterogeneous cube with alternating positive and
negative real permittivity and strictly positive loss. A constant complex
field and corresponding cellwise source provide an analytical absorbed
power. Doubling the load without changing the field must break source
power closure while preserving operator/cell agreement. A second test
covers zero power and invalid nonzero PEC coefficients.

This is a discrete energy diagnostic, not a scattering energy-balance
benchmark, a port/PML power audit, or proof of continuum accuracy. The
caller must provide the same mesh, material distribution, frequency and
homogeneous boundary assumptions used to assemble the system.
