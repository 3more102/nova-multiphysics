#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <span>
#include <vector>

namespace nova::fem::detail {

// One constant null mode per connected, unanchored continuous scalar H1 space.
class DiffusionNullspace {
public:
    DiffusionNullspace(std::size_t dof_count, std::size_t element_count,
        const std::function<std::span<const std::size_t>(std::size_t)>& element_dofs);

    void add_basis_integral(std::size_t dof, long double weight) {
        basis_integrals_[dof] += weight;
    }
    void record_load(std::size_t dof, double load) {
        components_[component_of_[dof]].absolute_load += std::abs(static_cast<long double>(load));
    }
    void anchor(std::size_t dof) { components_[component_of_[dof]].anchored = true; }

    void impose_gauges(std::span<const double> rhs, std::vector<std::optional<double>>& prescribed);
    [[nodiscard]] std::vector<double> physical_values(const std::vector<double>& iterate) const;

private:
    struct Component {
        std::vector<std::size_t> dofs;
        long double absolute_load{};
        bool anchored{};
    };
    std::vector<Component> components_;
    std::vector<std::size_t> component_of_;
    std::vector<long double> basis_integrals_;
};

} // namespace nova::fem::detail
