#include "heat_time_integration.hpp"

#include "nova/linalg/conjugate_gradient.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace nova::fem::detail {
namespace {

struct ConstrainedSystem {
    linalg::CsrMatrix matrix;
    std::vector<linalg::Triplet> couplings;
};

[[nodiscard]] ConstrainedSystem constrain(
    const linalg::CsrMatrix& full, const std::vector<bool>& constrained
) {
    if (std::none_of(constrained.begin(), constrained.end(), [](bool flag) { return flag; })) {
        return {full, {}};
    }
    std::vector<linalg::Triplet> entries, couplings;
    entries.reserve(full.non_zeros());
    full.for_each_nonzero([&](std::size_t i, std::size_t j, double value) {
        if (!constrained[i]) {
            if (constrained[j]) {
                couplings.push_back({i, j, value});
            } else {
                entries.push_back({i, j, value});
            }
        }
    });
    for (std::size_t i = 0; i < constrained.size(); ++i) {
        if (constrained[i]) {
            entries.push_back({i, i, 1.0});
        }
    }
    return {linalg::CsrMatrix::from_triplets(full.rows(), full.columns(), entries),
            std::move(couplings)};
}

struct HeatRate {
    std::vector<double> net;
    double source{};
    double neumann{};
    double robin{};
};

[[nodiscard]] HeatRate evaluate_rate(
    const linalg::CsrMatrix& stiffness, const HeatLoads& loads, const std::vector<double>& state
) {
    const auto conduction = stiffness.multiply(state);
    const auto robin = loads.robin_matrix.multiply(state);
    HeatRate rate;
    rate.net.resize(state.size());
    for (std::size_t i = 0; i < state.size(); ++i) {
        const double exchange = loads.robin_ambient_rhs[i] - robin[i];
        rate.source += loads.volume_rhs[i];
        rate.neumann += loads.neumann_rhs[i];
        rate.robin += exchange;
        rate.net[i] = loads.volume_rhs[i] + loads.neumann_rhs[i] + exchange - conduction[i];
    }
    return rate;
}

} // namespace

void validate_heat_time_grid(const TransientHeatTimeOptions& options) {
    const double dt = options.time_step;
    const double end = options.start_time + static_cast<double>(options.steps) * dt;
    if (options.steps == 0 || !std::isfinite(dt) || dt <= 0.0 ||
        !std::isfinite(options.start_time) || !std::isfinite(end) || end <= options.start_time) {
        throw std::invalid_argument("Transient time grid must be finite and advance.");
    }
    if (options.scheme != TransientTimeScheme::BackwardEuler &&
        options.scheme != TransientTimeScheme::CrankNicolson) {
        throw std::invalid_argument("Unsupported transient time integration scheme.");
    }
}

TransientHeatResult integrate_heat(
    const ScalarH1Operators& operators, const TransientHeatTimeOptions& options,
    std::vector<double> state, const std::vector<bool>& constrained,
    const std::function<std::vector<double>(double)>& dirichlet_at,
    const std::function<HeatLoads(double)>& loads_at,
    linalg::CachedConjugateGradient* solver_cache
) {
    const double dt = options.time_step;
    if (state.size() != constrained.size() ||
        !std::all_of(state.begin(), state.end(), [](double x) { return std::isfinite(x); })) {
        throw std::invalid_argument("Initial temperature must be finite at every DOF.");
    }
    const std::size_t n = state.size();
    const auto initial_values = dirichlet_at(options.start_time);
    for (std::size_t i = 0; i < n; ++i) {
        if (constrained[i]) {
            state[i] = initial_values[i];
        }
    }
    const bool trapezoidal = options.scheme == TransientTimeScheme::CrankNicolson;
    const double implicit_dt = trapezoidal ? dt / 2.0 : dt;
    const double explicit_dt = trapezoidal ? dt / 2.0 : 0.0;
    const auto base = linalg::CsrMatrix::linear_combination(
        operators.mass, 1.0, operators.stiffness, implicit_dt);
    std::optional<ConstrainedSystem> base_constrained;
    linalg::CachedConjugateGradient local_cg_cache;
    auto& cached_cg = solver_cache != nullptr ? *solver_cache : local_cg_cache;
    HeatRate previous;
    if (trapezoidal) {
        previous = evaluate_rate(operators.stiffness, loads_at(options.start_time), state);
    }

    auto mass_state = operators.mass.multiply(state);
    double heat = std::accumulate(mass_state.begin(), mass_state.end(), 0.0);
    TransientHeatResult result;
    result.final_time = options.start_time;
    result.linear_solvers.reserve(options.steps);
    result.diagnostics.reserve(options.steps);
    for (std::size_t step = 0; step < options.steps; ++step) {
        const double time = options.start_time + static_cast<double>(step + 1) * dt;
        if (time <= result.final_time || !std::isfinite(time)) {
            throw std::invalid_argument("A transient time step does not advance at this precision.");
        }
        const auto loads = loads_at(time);
        const auto prescribed = dirichlet_at(time);
        std::vector<double> rhs = mass_state;
        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] += implicit_dt * (loads.volume_rhs[i] + loads.neumann_rhs[i] +
                                     loads.robin_ambient_rhs[i]);
            if (trapezoidal) {
                rhs[i] += explicit_dt * previous.net[i];
            }
            if (!std::isfinite(rhs[i])) {
                throw std::overflow_error("Transient right-hand side overflowed.");
            }
        }
        // Reuse the constrained system on every step without a Robin matrix.
        std::optional<ConstrainedSystem> robin_system;
        if (loads.robin_matrix.non_zeros() != 0) {
            robin_system = constrain(linalg::CsrMatrix::linear_combination(
                base, 1.0, loads.robin_matrix, implicit_dt), constrained);
        } else if (!base_constrained) {
            base_constrained = constrain(base, constrained);
        }
        const auto& system = robin_system ? *robin_system : *base_constrained;
        for (const auto& entry : system.couplings) {
            rhs[entry.row] -= entry.value * prescribed[entry.column];
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (constrained[i]) {
                rhs[i] = prescribed[i];
            }
        }
        // Start identity rows at their prescribed values. Their residual and all
        // Krylov updates are then exactly zero, including across a restart.
        auto solver = cached_cg.solve(
            system.matrix, rhs, options.solver_options, prescribed);
        if (!solver.converged) {
            throw std::runtime_error("Transient linear solve failed at step " +
                                     std::to_string(step + 1) + ": " + solver.message);
        }

        auto next_mass = operators.mass.multiply(solver.solution);
        auto current = evaluate_rate(operators.stiffness, loads, solver.solution);
        TransientHeatStepDiagnostics diagnostic;
        diagnostic.time = time;
        diagnostic.heat_content = std::accumulate(next_mass.begin(), next_mass.end(), 0.0);
        diagnostic.source_energy = implicit_dt * current.source + explicit_dt * previous.source;
        diagnostic.neumann_energy = implicit_dt * current.neumann + explicit_dt * previous.neumann;
        diagnostic.robin_energy = implicit_dt * current.robin + explicit_dt * previous.robin;
        for (std::size_t i = 0; i < n; ++i) {
            if (constrained[i]) {
                double reaction = next_mass[i] - mass_state[i] - implicit_dt * current.net[i];
                if (trapezoidal) {
                    reaction -= explicit_dt * previous.net[i];
                }
                diagnostic.dirichlet_energy += reaction;
            }
        }
        diagnostic.balance_residual = diagnostic.heat_content - heat - diagnostic.source_energy
            - diagnostic.neumann_energy - diagnostic.robin_energy - diagnostic.dirichlet_energy;
        if (!std::isfinite(diagnostic.balance_residual)) {
            throw std::overflow_error("Transient heat budget is non-finite.");
        }
        heat = diagnostic.heat_content;
        mass_state = std::move(next_mass);
        if (trapezoidal) {
            previous = std::move(current);
        }
        result.diagnostics.push_back(diagnostic);
        result.linear_solvers.push_back(std::move(solver));
        result.final_time = time;
    }
    result.temperature = result.linear_solvers.back().solution;
    return result;
}

} // namespace nova::fem::detail
