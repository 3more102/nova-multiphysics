# Uniform Bragg grating reference

The header-only C++20 API `nova/physics/optics_bragg.hpp` solves the
lossless uniform two-mode coupled-envelope boundary problem

    A' = i delta A + i kappa B
    B' = -i kappa A - i delta B
    A(0) = 1, B(L) = 0

where detuning and real coupling are in inverse metres, and L is in metres.
The amplitudes are r=B(0), t=A(L), excluding the carrier phase.
For dimensionless d=delta L and k=kappa L, the matrix exponential gives
r=i k S/(C-i d S), t=1/(C-i d S), where C=cosh(x), S=sinh(x)/x
and x=sqrt(k^2-d^2) in the stopband. In the passband use cos and sinc
with x=sqrt(d^2-k^2); at the edge C=S=1. In the stopband, factor exp(x)
from the denominator and numerator to avoid growing exponentials.
R=|r|^2 and T=|t|^2 conserve unit incident power.

Regression cases cover the exact resonance R=tanh(k)^2, band edge,
passband null, strong grating, zero length/coupling, symmetry and
401 sampled detunings. The standalone test compiled and passed with GCC
C++20 and -Wall -Wextra -Wpedantic -Werror. Repository-wide CI is pending.

Scope: analytical coupled-mode reference only. Not a full-wave FEM solver,
not a model of arbitrary gratings, chirp, apodization, loss, polarization
conversion, radiation or dispersion. Coupling and detuning must be
supplied by the caller. Extreme phase values may be numerically unresolved.
