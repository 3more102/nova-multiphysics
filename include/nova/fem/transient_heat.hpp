#pragma once

#include "nova/linalg/csr_matrix.hpp"
#include "nova/linalg/solver_types.hpp"

#include <cstddef>
#include <vector>

namespace nova::fem {

enum class TransientTimeScheme { BackwardEuler, CrankNicolson };

struct ScalarH1Operators {
    linalg::CsrMatrix mass;
    linalg::CsrMatrix stiffness;
};

struct TransientHeatTimeOptions {
    TransientTimeScheme scheme{TransientTimeScheme::BackwardEuler};
    double time_step{0.1};
    double start_time{0.0};
    std::size_t steps{1U};
    linalg::SolverOptions solver_options{};
};

// Each energy term is integrated over the step, positive into the domain.
// The Dirichlet term is the sum of constrained-row reactions.
struct TransientHeatStepDiagnostics {
    double time{};
    double heat_content{};
    double source_energy{};
    double neumann_energy{};
    double robin_energy{};
    double dirichlet_energy{};
    double balance_residual{};
};

struct TransientHeatResult {
    std::vector<double> temperature;
    double final_time{};
    std::vector<linalg::SolverResult> linear_solvers;
    std::vector<TransientHeatStepDiagnostics> diagnostics;
};

} // namespace nova::fem
