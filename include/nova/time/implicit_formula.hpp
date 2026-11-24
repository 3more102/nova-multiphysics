#pragma once

namespace nova::time {
enum class ImplicitScheme { BackwardEuler, CrankNicolson, BDF2 };

// dH/dt = a*(H_new-H_old) - b*(H_old-H_older).
// theta weights the new spatial residual. BDF2 starts with trapezoidal integration.
struct ImplicitFormula {
    double a{}, b{}, theta{};
};
[[nodiscard]] ImplicitFormula implicit_formula(ImplicitScheme scheme, double dt, double previous_dt = 0);
}
