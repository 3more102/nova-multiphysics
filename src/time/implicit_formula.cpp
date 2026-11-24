#include "nova/time/implicit_formula.hpp"
#include <cmath>
#include <stdexcept>

namespace nova::time {
ImplicitFormula implicit_formula(ImplicitScheme scheme, double dt, double previous_dt) {
    if (!std::isfinite(dt) || !(dt > 0) || !std::isfinite(1/dt) ||
        !std::isfinite(previous_dt) || previous_dt < 0)
        throw std::invalid_argument("Invalid implicit time increment.");
    switch (scheme) {
    case ImplicitScheme::BackwardEuler: return {1/dt, 0, 1};
    case ImplicitScheme::CrankNicolson: return {1/dt, 0, 0.5};
    case ImplicitScheme::BDF2:
        if (previous_dt == 0) return {1/dt, 0, 0.5};
        const double fraction = dt >= previous_dt ? 1/(1+previous_dt/dt) :
            (dt/previous_dt)/(1+dt/previous_dt);
        const ImplicitFormula result{(1+fraction)/dt, fraction/previous_dt, 1};
        if (!std::isfinite(result.a) || !std::isfinite(result.b))
            throw std::overflow_error("Variable-step BDF2 coefficients overflow.");
        return result;
    }
    throw std::invalid_argument("Unknown time integration scheme.");
}
}
