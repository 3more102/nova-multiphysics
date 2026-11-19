#pragma once

#include "nova/fem/transient_heat.hpp"
#include "nova/linalg/conjugate_gradient.hpp"

#include <functional>

namespace nova::fem::detail {

struct HeatLoads {
    linalg::CsrMatrix robin_matrix;
    std::vector<double> neumann_rhs;
    std::vector<double> robin_ambient_rhs;
    std::vector<double> volume_rhs;
};

void validate_heat_time_grid(const TransientHeatTimeOptions& options);

[[nodiscard]] TransientHeatResult integrate_heat(
    const ScalarH1Operators& operators, const TransientHeatTimeOptions& options,
    std::vector<double> initial_state, const std::vector<bool>& constrained,
    const std::function<std::vector<double>(double)>& dirichlet_at,
    const std::function<HeatLoads(double)>& loads_at,
    linalg::CachedConjugateGradient* solver_cache = nullptr
);

} // namespace nova::fem::detail
