#include "nova/physics/waveguide_modes_2d.hpp"

#include "nova/linalg/csr_matrix.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::physics {
namespace {

using EdgeKey = std::pair<std::size_t, std::size_t>;

[[nodiscard]] std::vector<bool> exterior_vertex_mask(const fem::TriangleMesh2D& mesh) {
    std::map<EdgeKey, std::size_t> incidence;
    for (const auto& element : mesh.elements) {
        for (std::size_t local = 0; local < 3; ++local) {
            auto first = element[local];
            auto second = element[(local + 1U) % 3U];
            if (first > second) std::swap(first, second);
            ++incidence[{first, second}];
        }
    }

    std::vector<bool> boundary(mesh.nodes.size(), false);
    for (const auto& [edge, count] : incidence) {
        if (count == 1U) {
            boundary[edge.first] = true;
            boundary[edge.second] = true;
        } else if (count != 2U) {
            throw std::invalid_argument("Invalid triangle-edge incidence in waveguide mesh.");
        }
    }
    return boundary;
}

[[nodiscard]] bool positive_finite(const double value) {
    return std::isfinite(value) && value > 0.0;
}

} // namespace

WaveguideModesResult2D solve_waveguide_modes_2d(
    const WaveguideModeProblem2D& problem,
    const std::size_t max_modes,
    const linalg::SymmetricEigenOptions& eigen_options) {
    if (max_modes == 0U) {
        throw std::invalid_argument("At least one waveguide mode must be requested.");
    }
    problem.mesh.validate();
    if (problem.materials.size() != problem.mesh.elements.size()) {
        throw std::invalid_argument("Waveguide material count must match triangle count.");
    }
    if (!positive_finite(problem.vacuum_wavelength_m) ||
        !positive_finite(problem.cladding_index)) {
        throw std::invalid_argument("Wavelength and cladding index must be finite and positive.");
    }
    if (problem.polarization != WaveguidePolarization2D::TE_Ez &&
        problem.polarization != WaveguidePolarization2D::TM_Hz) {
        throw std::invalid_argument("Unsupported scalar waveguide polarization.");
    }

    constexpr double pi = std::numbers::pi_v<double>;
    const double k0 = 2.0 * pi / problem.vacuum_wavelength_m;
    const double k0_squared = k0 * k0;
    if (!positive_finite(k0) || !positive_finite(k0_squared)) {
        throw std::overflow_error("Vacuum wavenumber exceeds the supported numeric range.");
    }

    double maximum_material_index = 0.0;
    for (const auto& material : problem.materials) {
        if (!positive_finite(material.relative_permittivity) ||
            !positive_finite(material.relative_permeability)) {
            throw std::invalid_argument(
                "Waveguide materials must have finite positive real permittivity and permeability.");
        }
        const double index = std::sqrt(material.relative_permittivity) *
                             std::sqrt(material.relative_permeability);
        if (!positive_finite(index)) {
            throw std::overflow_error("Waveguide material index exceeds the supported numeric range.");
        }
        maximum_material_index = std::max(maximum_material_index, index);
    }
    const double index_scale = std::max({1.0, problem.cladding_index, maximum_material_index});
    const double index_tolerance =
        64.0 * std::numeric_limits<double>::epsilon() * index_scale;
    if (problem.cladding_index > maximum_material_index + index_tolerance) {
        throw std::invalid_argument("Cladding index cannot exceed the maximum material index.");
    }

    const auto boundary = exterior_vertex_mask(problem.mesh);
    std::vector<bool> constrained(problem.mesh.nodes.size(), false);
    if (problem.dirichlet_boundary_nodes.empty()) {
        constrained = boundary;
    } else {
        for (const std::size_t node : problem.dirichlet_boundary_nodes) {
            if (node >= problem.mesh.nodes.size()) {
                throw std::out_of_range("Waveguide Dirichlet node is outside the mesh.");
            }
            if (!boundary[node]) {
                throw std::invalid_argument("Waveguide Dirichlet constraints must lie on the exterior boundary.");
            }
            if (constrained[node]) {
                throw std::invalid_argument("Duplicate waveguide Dirichlet node.");
            }
            constrained[node] = true;
        }
    }

    constexpr std::size_t invalid_index = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> free_index(problem.mesh.nodes.size(), invalid_index);
    std::size_t free_count = 0;
    WaveguideModesResult2D result;
    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        if (constrained[node]) {
            ++result.constrained_nodes;
        } else {
            free_index[node] = free_count++;
        }
    }
    result.free_nodes = free_count;
    if (free_count == 0U) {
        throw std::invalid_argument("Waveguide constraints leave no free scalar degrees of freedom.");
    }

    // No material contrast above the declared cladding means no bound modes.
    if (maximum_material_index <= problem.cladding_index + index_tolerance) {
        result.converged = true;
        return result;
    }

    if (problem.mesh.elements.size() >
        std::numeric_limits<std::size_t>::max() / 9U) {
        throw std::length_error("Waveguide element count exceeds sparse assembly limits.");
    }
    std::vector<linalg::Triplet> operator_entries;
    std::vector<linalg::Triplet> mass_entries;
    std::vector<double> element_areas;
    std::vector<double> element_mass_coefficients;
    operator_entries.reserve(problem.mesh.elements.size() * 9U);
    mass_entries.reserve(problem.mesh.elements.size() * 9U);
    element_areas.reserve(problem.mesh.elements.size());
    element_mass_coefficients.reserve(problem.mesh.elements.size());

    for (std::size_t element_index = 0;
         element_index < problem.mesh.elements.size(); ++element_index) {
        const auto& element = problem.mesh.elements[element_index];
        const auto& p0 = problem.mesh.nodes[element[0]];
        const auto& p1 = problem.mesh.nodes[element[1]];
        const auto& p2 = problem.mesh.nodes[element[2]];
        const double determinant =
            (p1.x - p0.x) * (p2.y - p0.y) -
            (p2.x - p0.x) * (p1.y - p0.y);
        const double area = 0.5 * std::abs(determinant);
        if (!positive_finite(area) || !std::isfinite(determinant)) {
            throw std::invalid_argument("Invalid triangle geometry in waveguide assembly.");
        }

        const std::array<std::array<double, 2>, 3> gradients{{
            {{(p1.y - p2.y) / determinant, (p2.x - p1.x) / determinant}},
            {{(p2.y - p0.y) / determinant, (p0.x - p2.x) / determinant}},
            {{(p0.y - p1.y) / determinant, (p1.x - p0.x) / determinant}}
        }};

        const auto& material = problem.materials[element_index];
        double gradient_coefficient = 0.0;
        double potential_coefficient = 0.0;
        double mass_coefficient = 0.0;
        if (problem.polarization == WaveguidePolarization2D::TE_Ez) {
            gradient_coefficient = 1.0 / material.relative_permeability;
            potential_coefficient = material.relative_permittivity;
            mass_coefficient = 1.0 / material.relative_permeability;
        } else {
            gradient_coefficient = 1.0 / material.relative_permittivity;
            potential_coefficient = material.relative_permeability;
            mass_coefficient = 1.0 / material.relative_permittivity;
        }

        if (!positive_finite(gradient_coefficient) ||
            !positive_finite(potential_coefficient) ||
            !positive_finite(mass_coefficient)) {
            throw std::overflow_error("Waveguide material weighting exceeds the supported numeric range.");
        }
        const double mass_scale = area / 12.0;
        const double potential_scale = k0_squared * potential_coefficient * mass_scale;
        const double gradient_scale = area * gradient_coefficient;
        const double weighted_mass_scale = mass_coefficient * mass_scale;
        if (!positive_finite(mass_scale) || !positive_finite(potential_scale) ||
            !positive_finite(gradient_scale) || !positive_finite(weighted_mass_scale)) {
            throw std::overflow_error("Waveguide element coefficient overflow.");
        }
        element_areas.push_back(area);
        element_mass_coefficients.push_back(mass_coefficient);

        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                const double gradient_dot =
                    gradients[i][0] * gradients[j][0] +
                    gradients[i][1] * gradients[j][1];
                const double shape_mass = (i == j ? 2.0 : 1.0);
                const double local_operator =
                    potential_scale * shape_mass - gradient_scale * gradient_dot;
                const double local_mass = weighted_mass_scale * shape_mass;
                if (!std::isfinite(local_operator) || !std::isfinite(local_mass)) {
                    throw std::overflow_error("Waveguide local matrix entry overflow.");
                }

                const std::size_t row = free_index[element[i]];
                const std::size_t column = free_index[element[j]];
                if (row != invalid_index && column != invalid_index) {
                    operator_entries.push_back({row, column, local_operator});
                    mass_entries.push_back({row, column, local_mass});
                }
            }
        }
    }

    const auto wave_operator = linalg::CsrMatrix::from_triplets(
        free_count, free_count, operator_entries);
    const auto mass = linalg::CsrMatrix::from_triplets(
        free_count, free_count, mass_entries);
    const auto eigen = linalg::solve_symmetric_generalized_eigen(
        wave_operator, mass, eigen_options);

    result.jacobi_sweeps = eigen.sweeps;
    result.converged = eigen.converged;

    const double cladding_beta = k0 * problem.cladding_index;
    const double maximum_beta = k0 * maximum_material_index;
    const double cutoff_squared = cladding_beta * cladding_beta;
    const double maximum_squared = maximum_beta * maximum_beta;
    if (!std::isfinite(cutoff_squared) || !std::isfinite(maximum_squared)) {
        throw std::overflow_error("Waveguide propagation-constant bound overflow.");
    }
    const double spectral_scale = std::max(1.0, maximum_squared);
    const double spectral_tolerance = std::max(
        128.0 * std::numeric_limits<double>::epsilon() * spectral_scale,
        100.0 * eigen_options.relative_tolerance * spectral_scale);

    for (auto it = eigen.eigenpairs.rbegin();
         it != eigen.eigenpairs.rend() && result.modes.size() < max_modes; ++it) {
        const double beta_squared = it->eigenvalue;
        if (!std::isfinite(beta_squared)) {
            throw std::runtime_error("Waveguide eigensolver returned a nonfinite eigenvalue.");
        }
        if (beta_squared <= cutoff_squared + spectral_tolerance) {
            continue; // radiation continuum and cutoff modes are not bound modes
        }
        if (beta_squared > maximum_squared + spectral_tolerance) {
            throw std::runtime_error("Waveguide eigenvalue exceeds the maximum material light line.");
        }
        if (beta_squared >= maximum_squared - spectral_tolerance) {
            continue; // exclude the strict core-index endpoint
        }

        const double beta = std::sqrt(beta_squared);
        const double effective_index = beta / k0;
        if (!positive_finite(beta) ||
            !(effective_index > problem.cladding_index &&
              effective_index < maximum_material_index)) {
            continue;
        }

        WaveguideMode2D mode;
        mode.effective_index = effective_index;
        mode.propagation_constant_per_m = beta;
        mode.propagation_constant_squared_per_m2 = beta_squared;
        mode.relative_residual = it->relative_residual;
        mode.scalar_field.assign(problem.mesh.nodes.size(), 0.0);
        for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
            if (free_index[node] != invalid_index) {
                mode.scalar_field[node] = it->eigenvector[free_index[node]];
            }
        }

        const auto largest = std::max_element(
            mode.scalar_field.begin(), mode.scalar_field.end(),
            [](const double lhs, const double rhs) {
                return std::abs(lhs) < std::abs(rhs);
            });
        if (largest != mode.scalar_field.end() && *largest < 0.0) {
            for (double& value : mode.scalar_field) value = -value;
        }

        double total_weighted_norm = 0.0;
        double core_weighted_norm = 0.0;
        for (std::size_t element_index = 0;
             element_index < problem.mesh.elements.size(); ++element_index) {
            const auto& element = problem.mesh.elements[element_index];
            const double u0 = mode.scalar_field[element[0]];
            const double u1 = mode.scalar_field[element[1]];
            const double u2 = mode.scalar_field[element[2]];
            const double field_sum = u0 + u1 + u2;
            const double quadratic_sum =
                u0 * u0 + u1 * u1 + u2 * u2 + field_sum * field_sum;
            const double weighted_integral =
                element_mass_coefficients[element_index] *
                element_areas[element_index] / 12.0 * quadratic_sum;
            if (!std::isfinite(weighted_integral) || weighted_integral < 0.0) {
                throw std::runtime_error("Invalid waveguide modal mass integral.");
            }
            total_weighted_norm += weighted_integral;
            const auto& material = problem.materials[element_index];
            const double material_index =
                std::sqrt(material.relative_permittivity) *
                std::sqrt(material.relative_permeability);
            if (material_index > problem.cladding_index + index_tolerance) {
                core_weighted_norm += weighted_integral;
            }
        }
        if (!positive_finite(total_weighted_norm) ||
            !std::isfinite(core_weighted_norm)) {
            throw std::runtime_error("Invalid waveguide modal normalization integral.");
        }
        mode.core_confinement_factor = std::clamp(
            core_weighted_norm / total_weighted_norm, 0.0, 1.0);
        result.modes.push_back(std::move(mode));
    }
    return result;
}

} // namespace nova::physics
