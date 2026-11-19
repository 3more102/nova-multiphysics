#include "constrained_diffusion.hpp"

#include "nova/linalg/conjugate_gradient.hpp"

#include <cmath>
#include <stdexcept>

namespace nova::fem::detail {

linalg::SolverResult solve_constrained_diffusion(
    std::vector<linalg::Triplet> entries, std::vector<double> rhs,
    const std::vector<std::optional<double>>& prescribed, const linalg::SolverOptions& options
) {
    // Compact the free/free block in the original storage as we eliminate
    // prescribed columns. This avoids a second element-sized triplet list.
    std::size_t kept = 0;
    for (const auto& entry : entries) {
        if (!std::isfinite(entry.value)) {
            throw std::overflow_error("Non-finite diffusion entry before constraint elimination.");
        }
        if (prescribed[entry.row]) {
            continue;
        }
        if (prescribed[entry.column]) {
            rhs[entry.row] -= entry.value * *prescribed[entry.column];
            if (!std::isfinite(rhs[entry.row])) {
                throw std::overflow_error("Diffusion Dirichlet RHS correction overflowed.");
            }
        } else {
            entries[kept++] = entry;
        }
    }
    entries.resize(kept);
    std::vector<double> initial(rhs.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        if (prescribed[i]) {
            entries.push_back({i, i, 1.0});
            initial[i] = rhs[i] = *prescribed[i];
        }
    }
    const auto matrix = linalg::CsrMatrix::from_triplets(rhs.size(), rhs.size(), entries);
    return linalg::solve_conjugate_gradient(matrix, rhs, options, initial);
}

} // namespace nova::fem::detail
