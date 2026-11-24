# Discrete local optical power audit

The header `nova/physics/optical_helmholtz_power_audit.hpp` provides an experimental postprocessor for the normal-incidence P1 Helmholtz solver. It reports normalized signed rightward flux at both endpoints of every element, local absorption residual, maximum interior interface flux jump, and exterior endpoint residuals. Supply exactly the mesh, permittivity, incident index and wavelength used for the field solve.

Let h=k0*(x_right-x_left), a=E_left, b=E_right, and e=relative permittivity. The element matrix has diagonal 1/h-e*h/3 and off-diagonal -1/h-e*h/6. Its endpoint reactions fL and fR give normalized powers -Im(conj(a)*fL)/n_inc and Im(conj(b)*fR)/n_inc respectively. These are weak-form discrete powers, not samples of the raw P1 derivative. Their difference equals Im(e)*h*(|a|²+|b|²+Re(conj(a)*b))/(3*n_inc). The audit compares this difference with the stored cell absorptance instead of defining flux from that value.

At a converged discrete solution, shared-node reactions cancel and the exterior powers equal 1-R and T. A small audit residual is an algebraic consistency check, not evidence of mesh convergence or physical accuracy. No automatic spatial refinement, oblique TE/TM support or 3D Maxwell coupling is provided here.

The CTest target `nova_helmholtz_power_audit_tests` checks a nonuniform mesh with transparent, absorbing and metal-like regions and unequal exterior indices. It also deliberately modifies a nodal field and cell loss to check detection of inconsistent data.
