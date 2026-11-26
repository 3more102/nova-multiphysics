#include "nova/physics/floating_electrodes_3d.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
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


FloatingElectrodePotentials3D
solve_affine_floating_electrode_potentials_3d(
    const CapacitanceMatrix3D& capacitance,
    const std::span<const double> zero_voltage_reaction_offsets_coulombs,
    const std::span<const FixedElectrodePotential3D> fixed_electrodes,
    const std::span<const FloatingElectrodeCharge3D> floating_electrodes
) {
    const std::size_t count = capacitance.electrode_ids.size();
    if (zero_voltage_reaction_offsets_coulombs.size() != count) {
        throw std::invalid_argument(
            "Affine floating-electrode offsets must match the terminal count.");
    }
    for (const double charge : zero_voltage_reaction_offsets_coulombs) {
        if (!std::isfinite(charge)) {
            throw std::invalid_argument(
                "Affine floating-electrode offsets must be finite.");
        }
    }

    // The homogeneous terminal matrix is unchanged by impressed sources.
    // Solve C*V = Q_requested - Q_zero, using the existing pivot validation.
    std::vector<FloatingElectrodeCharge3D> shifted;
    shifted.reserve(floating_electrodes.size());
    for (const auto& terminal : floating_electrodes) {
        if (!std::isfinite(terminal.charge_coulombs)) {
            throw std::invalid_argument(
                "Affine floating-electrode prescribed charges must be finite.");
        }
        const auto found = std::lower_bound(
            capacitance.electrode_ids.begin(), capacitance.electrode_ids.end(),
            terminal.electrode_id);
        if (found == capacitance.electrode_ids.end() ||
            *found != terminal.electrode_id) {
            throw std::invalid_argument(
                "Unknown electrode ID in affine floating-electrode inputs.");
        }
        const auto index = static_cast<std::size_t>(
            found - capacitance.electrode_ids.begin());
        shifted.push_back({
            terminal.electrode_id,
            finite_double(static_cast<long double>(terminal.charge_coulombs) -
                          zero_voltage_reaction_offsets_coulombs[index])
        });
    }

    auto result = solve_floating_electrode_potentials_3d(
        capacitance, fixed_electrodes, shifted);
    for (std::size_t i = 0; i < count; ++i) {
        result.estimated_charge_coulombs[i] = finite_double(
            static_cast<long double>(result.estimated_charge_coulombs[i]) +
            zero_voltage_reaction_offsets_coulombs[i]);
    }

    // Residuals are defined against the *actual* requested total charge,
    // not the intermediate source-subtracted terminal right-hand side.
    result.max_floating_charge_residual_coulombs = 0.0;
    for (const auto& terminal : floating_electrodes) {
        const auto index = static_cast<std::size_t>(
            std::lower_bound(result.electrode_ids.begin(),
                             result.electrode_ids.end(),
                             terminal.electrode_id) -
            result.electrode_ids.begin());
        result.max_floating_charge_residual_coulombs = std::max(
            result.max_floating_charge_residual_coulombs,
            finite_double(std::abs(
                static_cast<long double>(result.estimated_charge_coulombs[index]) -
                static_cast<long double>(terminal.charge_coulombs))));
    }
    return result;
}


ReconstructedFloatingElectrostaticField3D
reconstruct_floating_electrostatic_field_3d(
    const Electrostatics3DProblem& dielectric_problem,
    const std::span<const ElectrodeBoundaryGroup3D> electrode_groups,
    const std::span<const FixedElectrodePotential3D> fixed_electrodes,
    const std::span<const FloatingElectrodeCharge3D> floating_electrodes,
    const linalg::SolverOptions& solver_options
) {
    // Capacitance extraction intentionally removes all impressed sources.
    // A field reconstruction from that matrix would therefore be misleading
    // if the caller supplied charge or nonzero displacement loads.
    if (dielectric_problem.free_charge_density ||
        !dielectric_problem.cell_free_charge_density.empty() ||
        !dielectric_problem.displacement_flux_faces.empty()) {
        throw std::invalid_argument(
            "Charge-to-field reconstruction requires a source-free dielectric "
            "and no prescribed displacement boundary loads.");
    }

    ReconstructedFloatingElectrostaticField3D result;
    result.capacitance = extract_capacitance_matrix_3d(
        dielectric_problem, electrode_groups, solver_options);
    result.terminals = solve_floating_electrode_potentials_3d(
        result.capacitance, fixed_electrodes, floating_electrodes);

    // The extractor validates the entire group partition, including
    // distinct exterior faces and no cross-electrode shared vertices.
    // Canonical face keys handle any input face winding.
    std::set<fem::BoundaryFace3D> electrode_faces;
    std::set<std::size_t> electrode_nodes;
    std::vector<fem::BoundaryFace3D> all_faces;
    std::vector<int> all_ids;
    for (const auto& group : electrode_groups) {
        for (const auto& face : group.faces) {
            auto canonical = face;
            std::sort(canonical.begin(), canonical.end());
            electrode_faces.insert(canonical);
            for (const auto node : face) electrode_nodes.insert(node);
            all_faces.push_back(face);
            all_ids.push_back(group.electrode_id);
        }
    }

    Electrostatics3DProblem field_problem = dielectric_problem;
    field_problem.potential_constraints.clear();
    field_problem.potential_faces.clear();

    // Keep non-electrode reference locations, but ground them at exactly 0 V,
    // just as extract_capacitance_matrix_3d does for each basis solve.
    for (const auto& condition : dielectric_problem.potential_constraints) {
        if (!electrode_nodes.contains(condition.node)) {
            field_problem.potential_constraints.push_back({
                condition.node, 0.0
            });
        }
    }
    const fem::ScalarField3D zero = [](double, double, double) {
        return 0.0;
    };
    for (const auto& condition : dielectric_problem.potential_faces) {
        auto canonical = condition.face;
        std::sort(canonical.begin(), canonical.end());
        if (!electrode_faces.contains(canonical)) {
            field_problem.potential_faces.push_back({
                condition.face, zero
            });
        }
    }

    for (const auto& group : electrode_groups) {
        const auto it = std::lower_bound(
            result.terminals.electrode_ids.begin(),
            result.terminals.electrode_ids.end(), group.electrode_id);
        if (it == result.terminals.electrode_ids.end() ||
            *it != group.electrode_id) {
            throw std::logic_error("Missing inferred electrode potential.");
        }
        const auto index = static_cast<std::size_t>(
            it - result.terminals.electrode_ids.begin());
        const double volts = result.terminals.potential_volts.at(index);
        const auto constant_volts = [volts](double, double, double) {
            return volts;
        };
        for (const auto& face : group.faces) {
            field_problem.potential_faces.push_back({
                face, constant_volts
            });
        }
    }

    result.field = solve_electrostatics_3d(field_problem, solver_options);
    if (!result.field.linear_solver.converged) {
        throw std::runtime_error(
            "Floating-electrode reconstructed field did not converge.");
    }
    result.electrode_charges = evaluate_electrode_charges_3d(
        field_problem, result.field, all_faces, all_ids);
    if (result.electrode_charges.size() != result.terminals.electrode_ids.size()) {
        throw std::runtime_error("Reconstructed electrode-charge set incomplete.");
    }
    for (const auto& prescribed : floating_electrodes) {
        const auto found = std::lower_bound(
            result.electrode_charges.begin(), result.electrode_charges.end(),
            prescribed.electrode_id,
            [](const ElectrodeCharge3D& a, const int id) {
                return a.electrode_id < id;
            });
        if (found == result.electrode_charges.end() ||
            found->electrode_id != prescribed.electrode_id) {
            throw std::logic_error("Missing reconstructed floating-electrode charge.");
        }
        const auto difference = std::abs(
            static_cast<long double>(found->conductor_charge_coulombs) -
            static_cast<long double>(prescribed.charge_coulombs));
        result.max_floating_charge_mismatch_coulombs = std::max(
            result.max_floating_charge_mismatch_coulombs,
            finite_double(difference));
    }
    return result;
}

} // namespace nova::physics
