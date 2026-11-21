#include "nova/physics/floating_electrodes_3d.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nova::physics {
namespace {

[[nodiscard]] std::size_t electrode_index(
    const std::vector<int>& sorted_ids, const int id
) {
    const auto found = std::lower_bound(sorted_ids.begin(), sorted_ids.end(), id);
    if (found == sorted_ids.end() || *found != id) {
        throw std::invalid_argument("Unknown electrode ID in floating-electrode inputs.");
    }
    return static_cast<std::size_t>(found - sorted_ids.begin());
}

[[nodiscard]] double finite_double(const long double value) {
    if (!std::isfinite(value) ||
        std::abs(value) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Floating-electrode calculation exceeded finite double range.");
    }
    return static_cast<double>(value);
}

} // namespace

FloatingElectrodePotentials3D solve_floating_electrode_potentials_3d(
    const CapacitanceMatrix3D& capacitance,
    const std::span<const FixedElectrodePotential3D> fixed_electrodes,
    const std::span<const FloatingElectrodeCharge3D> floating_electrodes
) {
    const auto& ids = capacitance.electrode_ids;
    const auto& c = capacitance.matrix_farads;
    const std::size_t n = ids.size();
    if (n == 0U || n > c.max_size() / n || c.size() != n * n) {
        throw std::invalid_argument("Capacitance must be a nonempty square matrix.");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (ids[i] < 0 || (i > 0U && ids[i - 1U] >= ids[i])) {
            throw std::invalid_argument("Capacitance electrode IDs must be strictly increasing and nonnegative.");
        }
    }
    for (const double entry : c) {
        if (!std::isfinite(entry)) {
            throw std::invalid_argument("Capacitance matrix entries must be finite.");
        }
    }

    // 0 = unspecified, 1 = fixed voltage, 2 = floating with prescribed Q.
    std::vector<int> state(n, 0);
    std::vector<long double> v(n, 0.0L);
    std::vector<long double> prescribed_charge(n, 0.0L);
    std::vector<std::size_t> free_indices;
    free_indices.reserve(floating_electrodes.size());

    for (const auto& item : fixed_electrodes) {
        const std::size_t i = electrode_index(ids, item.electrode_id);
        if (state[i] != 0 || !std::isfinite(item.potential_volts)) {
            throw std::invalid_argument("Fixed electrode voltage is duplicate or non-finite.");
        }
        state[i] = 1;
        v[i] = static_cast<long double>(item.potential_volts);
    }
    for (const auto& item : floating_electrodes) {
        const std::size_t i = electrode_index(ids, item.electrode_id);
        if (state[i] != 0 || !std::isfinite(item.charge_coulombs)) {
            throw std::invalid_argument("Floating electrode charge is duplicate or non-finite.");
        }
        state[i] = 2;
        prescribed_charge[i] = static_cast<long double>(item.charge_coulombs);
        free_indices.push_back(i);
    }
    if (std::any_of(state.begin(), state.end(),
                    [](const int value) { return value == 0; })) {
        throw std::invalid_argument("Every capacitance electrode must be assigned exactly once.");
    }

    const std::size_t m = free_indices.size();
    if (m > 0U) {
        if (m == std::numeric_limits<std::size_t>::max() ||
            m > std::vector<long double>{}.max_size() / (m + 1U)) {
            throw std::length_error("Floating-electrode matrix dimensions overflow.");
        }
        const std::size_t width = m + 1U;
        std::vector<long double> a(m * width, 0.0L);
        long double matrix_scale = 0.0L;

        for (std::size_t r = 0; r < m; ++r) {
            const std::size_t i = free_indices[r];
            long double rhs = prescribed_charge[i];
            for (std::size_t j = 0; j < n; ++j) {
                if (state[j] == 1) {
                    rhs -= static_cast<long double>(c[i * n + j]) * v[j];
                }
            }
            a[r * width + m] = rhs;
            for (std::size_t col = 0; col < m; ++col) {
                const long double value =
                    static_cast<long double>(c[i * n + free_indices[col]]);
                a[r * width + col] = value;
                matrix_scale = std::max(matrix_scale, std::abs(value));
            }
        }
        if (!std::isfinite(matrix_scale) || matrix_scale == 0.0L) {
            throw std::invalid_argument("Floating-electrode capacitance block is singular.");
        }

        // Gaussian elimination with row pivoting. The relative threshold
        // deliberately rejects an ill-conditioned charge-to-voltage inverse.
        constexpr long double pivot_fraction = 1.0e-12L;
        for (std::size_t p = 0; p < m; ++p) {
            std::size_t pivot_row = p;
            long double pivot_abs = std::abs(a[p * width + p]);
            for (std::size_t r = p + 1U; r < m; ++r) {
                const long double candidate = std::abs(a[r * width + p]);
                if (candidate > pivot_abs) {
                    pivot_abs = candidate;
                    pivot_row = r;
                }
            }
            if (!(pivot_abs > pivot_fraction * matrix_scale)) {
                throw std::invalid_argument("Floating-electrode capacitance block is singular or ill-conditioned.");
            }
            if (pivot_row != p) {
                for (std::size_t col = 0; col < width; ++col) {
                    std::swap(a[p * width + col], a[pivot_row * width + col]);
                }
            }
            for (std::size_t r = p + 1U; r < m; ++r) {
                const long double factor = a[r * width + p] / a[p * width + p];
                a[r * width + p] = 0.0L;
                for (std::size_t col = p + 1U; col < width; ++col) {
                    a[r * width + col] -= factor * a[p * width + col];
                }
            }
        }
        for (std::size_t row = m; row-- > 0U;) {
            long double rhs = a[row * width + m];
            for (std::size_t col = row + 1U; col < m; ++col) {
                rhs -= a[row * width + col] * v[free_indices[col]];
            }
            v[free_indices[row]] = rhs / a[row * width + row];
            (void)finite_double(v[free_indices[row]]);
        }
    }

    FloatingElectrodePotentials3D result;
    result.electrode_ids = ids;
    result.potential_volts.reserve(n);
    result.estimated_charge_coulombs.reserve(n);
    for (const long double voltage : v) {
        result.potential_volts.push_back(finite_double(voltage));
    }
    for (std::size_t i = 0; i < n; ++i) {
        long double q = 0.0L;
        for (std::size_t j = 0; j < n; ++j) {
            q += static_cast<long double>(c[i * n + j]) * v[j];
        }
        result.estimated_charge_coulombs.push_back(finite_double(q));
        if (state[i] == 2) {
            result.max_floating_charge_residual_coulombs =
                std::max(result.max_floating_charge_residual_coulombs,
                         finite_double(std::abs(q - prescribed_charge[i])));
        }
    }
    return result;
}

} // namespace nova::physics
