#pragma once

#include "nova/fem/transient_heat_3d.hpp"
#include "nova/linalg/conjugate_gradient.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace nova::fem::detail {

// Per-integration prepared operators for the time-independent P1 tetrahedral
// heat model. Two bounded time-step slots let adaptive step-doubling retain
// the full-step and half-step systems independently.
class PreparedTransientHeat3D {
public:
    explicit PreparedTransientHeat3D(const TransientHeat3DProblem& problem);

    [[nodiscard]] TransientHeat3DResult solve(
        const std::vector<double>& initial_temperature,
        const TransientHeat3DOptions& options,
        std::size_t cache_slot = 0U
    );

    // Internal regression metrics; not part of the installed NOVA SDK.
    [[nodiscard]] std::size_t system_assemblies() const noexcept {
        return system_assemblies_;
    }
    [[nodiscard]] std::size_t factor_builds() const noexcept {
        return factor_builds_;
    }

private:
    struct StepSystemCache {
        bool initialized{false};
        double time_step{};
        TransientHeat3DScheme scheme{TransientHeat3DScheme::BackwardEuler};
        linalg::CsrMatrix lhs;
        linalg::CsrMatrix rhs_operator;
        linalg::CsrMatrix system;
        std::vector<double> boundary_shift;
        linalg::CachedConjugateGradient cached_cg;
    };

    std::size_t node_count_{};
    std::vector<std::optional<double>> prescribed_;
    linalg::CsrMatrix mass_;
    linalg::CsrMatrix stiffness_;
    linalg::CsrMatrix robin_;
    std::vector<double> volume_load_;
    std::vector<double> neumann_load_;
    std::vector<double> ambient_load_;
    std::vector<double> total_load_;
    std::vector<double> capacity_weights_;
    std::array<StepSystemCache, 2> step_systems_{};
    std::size_t system_assemblies_{};
    std::size_t factor_builds_{};
};

} // namespace nova::fem::detail
