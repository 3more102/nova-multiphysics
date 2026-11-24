#include "diffusion_nullspace.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace nova::fem::detail {
namespace {

class DisjointSet {
public:
    explicit DisjointSet(const std::size_t size)
        : parent_(size), rank_(size, 0U) {
        for (std::size_t node = 0; node < size; ++node) {
            parent_[node] = node;
        }
    }

    [[nodiscard]] std::size_t find(const std::size_t node) {
        if (parent_[node] != node) {
            parent_[node] = find(parent_[node]);
        }
        return parent_[node];
    }

    void unite(const std::size_t first, const std::size_t second) {
        std::size_t first_root = find(first);
        std::size_t second_root = find(second);
        if (first_root == second_root) {
            return;
        }
        if (rank_[first_root] < rank_[second_root]) {
            std::swap(first_root, second_root);
        }
        parent_[second_root] = first_root;
        if (rank_[first_root] == rank_[second_root]) {
            ++rank_[first_root];
        }
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<std::size_t> rank_;
};

} // namespace

DiffusionNullspace::DiffusionNullspace(
    std::size_t dof_count, std::size_t element_count,
    const std::function<std::span<const std::size_t>(std::size_t)>& element_dofs
) : component_of_(dof_count), basis_integrals_(dof_count, 0.0L) {
    DisjointSet sets(dof_count);
    for (std::size_t e = 0; e < element_count; ++e) {
        const auto local = element_dofs(e);
        for (const auto dof : local.subspan(1)) {
            sets.unite(local.front(), dof);
        }
    }
    // DOFs are numbered densely. The first DOF encountered gives each
    // component its stable identifier and its reproducible gauge vertex.
    std::vector<std::size_t> index_for_root(dof_count, dof_count);
    for (std::size_t dof = 0; dof < dof_count; ++dof) {
        auto& index = index_for_root[sets.find(dof)];
        if (index == dof_count) {
            index = components_.size();
            components_.emplace_back();
        }
        component_of_[dof] = index;
        components_[index].dofs.push_back(dof);
    }
}

void DiffusionNullspace::impose_gauges(
    std::span<const double> rhs, std::vector<std::optional<double>>& prescribed,
    bool require_anchors
) {
    for (std::size_t dof = 0; dof < prescribed.size(); ++dof) {
        if (prescribed[dof]) {
            anchor(dof);
        }
    }
    constexpr long double relative_tolerance = 1.0e-10L;
    for (const auto& c : components_) {
        if (c.anchored) {
            continue;
        }
        if (require_anchors) {
            throw std::invalid_argument("Absolute temperature requires a Dirichlet or positive Robin anchor "
                                        "on the component containing DOF " + std::to_string(c.dofs.front()));
        }
        long double net = 0.0L;
        for (const auto dof : c.dofs) {
            net += static_cast<long double>(rhs[dof]);
        }
        if (!std::isfinite(net) || !std::isfinite(c.absolute_load)) {
            throw std::overflow_error("Pure-Neumann component load accumulation overflowed.");
        }
        if (std::abs(net) > relative_tolerance * c.absolute_load) {
            throw std::invalid_argument("Incompatible pure-Neumann component containing DOF " +
                std::to_string(c.dofs.front()) + ": source and outward flux must balance.");
        }
        prescribed[c.dofs.front()] = 0.0;
    }
}

std::vector<double> DiffusionNullspace::physical_values(const std::vector<double>& iterate) const {
    auto values = iterate;
    for (const auto& c : components_) {
        if (c.anchored) {
            continue;
        }
        long double measure = 0.0L, integral = 0.0L;
        for (const auto dof : c.dofs) {
            measure += basis_integrals_[dof];
            integral += basis_integrals_[dof] * static_cast<long double>(values[dof]);
        }
        if (!(measure > 0.0L) || !std::isfinite(measure)) {
            throw std::overflow_error("Pure-Neumann component measure must be finite and positive.");
        }
        const double mean = static_cast<double>(integral / measure);
        if (!std::isfinite(mean)) {
            throw std::overflow_error("Pure-Neumann field mean overflowed.");
        }
        for (const auto dof : c.dofs) {
            values[dof] -= mean;
            if (!std::isfinite(values[dof])) {
                throw std::overflow_error("Pure-Neumann field normalization overflowed.");
            }
        }
    }
    return values;
}

} // namespace nova::fem::detail
