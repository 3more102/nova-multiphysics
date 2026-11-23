#include "nova/physics/electrothermal_3d.hpp"

#include "nova/fem/p1_flux_3d.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace nova::physics {
namespace {

void require_same_mesh(const fem::TetrahedralMesh3D& a,
                       const fem::TetrahedralMesh3D& b) {
    if (a.nodes.size() != b.nodes.size() ||
        a.elements.size() != b.elements.size()) {
        throw std::invalid_argument(
            "Electrothermal 3D requires identical electrical and thermal meshes.");
    }
    for (std::size_t i = 0; i < a.nodes.size(); ++i) {
        if (a.nodes[i].x != b.nodes[i].x ||
            a.nodes[i].y != b.nodes[i].y ||
            a.nodes[i].z != b.nodes[i].z) {
            throw std::invalid_argument(
                "Electrothermal 3D meshes have different coordinates.");
        }
    }
    for (std::size_t i = 0; i < a.elements.size(); ++i) {
        if (a.elements[i] != b.elements[i]) {
            throw std::invalid_argument(
                "Electrothermal 3D meshes have different element ordering.");
        }
    }
}

[[nodiscard]] double checked_double(const long double value, const char* message) {
    if (!std::isfinite(value) ||
        value > static_cast<long double>(std::numeric_limits<double>::max()) ||
        value < -static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error(message);
    }
    return static_cast<double>(value);
}

} // namespace

Electrothermal3DResult solve_one_way_electrothermal_3d(
    const Electrothermal3DProblem& problem,
    const linalg::SolverOptions& electrical_options,
    const linalg::SolverOptions& thermal_options
) {
    require_same_mesh(problem.electrical.mesh, problem.thermal.mesh);

    Electrothermal3DResult result;
    result.electrical_solution =
        fem::solve_scalar_diffusion_3d(problem.electrical, electrical_options);
    if (!result.electrical_solution.linear_solver.converged) {
        throw std::runtime_error("Electrothermal 3D electrical solve did not converge.");
    }
    // Flux averages use precisely the conductivity selection and quadrature
    // of the stationary P1 electrical stiffness assembly, scalar or SPD tensor.
    const auto fluxes = fem::evaluate_p1_tetrahedral_flux_3d(
        problem.electrical, result.electrical_solution.values);
    result.electric_field.reserve(fluxes.size());
    result.current_density.reserve(fluxes.size());
    result.joule_density.reserve(fluxes.size());

    long double integrated_power = 0.0L;
    for (std::size_t element = 0; element < fluxes.size(); ++element) {
        const auto& flux = fluxes[element];
        const fem::Point3D field{
            -flux.gradient.x, -flux.gradient.y, -flux.gradient.z
        };
        const long double x = static_cast<long double>(field.x) * flux.flux.x;
        const long double y = static_cast<long double>(field.y) * flux.flux.y;
        const long double z = static_cast<long double>(field.z) * flux.flux.z;
        const long double power_density = x + y + z;
        const long double scale = std::abs(x) + std::abs(y) + std::abs(z);
        // Positive-definite conductivity guarantees nonnegative dissipation.
        // Tolerate only last-bit cancellation for off-diagonal tensors.
        if (!std::isfinite(power_density) ||
            power_density < -256.0L *
                std::numeric_limits<long double>::epsilon() * scale) {
            throw std::overflow_error(
                "Electrothermal 3D Joule density is non-finite or negative.");
        }
        const long double nonnegative = std::max(0.0L, power_density);
        const double density = checked_double(
            nonnegative, "Electrothermal 3D Joule density exceeds double range.");
        const auto basis = fem::evaluate_tetrahedron_p1_physical_gradients_3d(
            problem.electrical.mesh, element);
        const long double volume =
            std::abs(static_cast<long double>(basis.signed_jacobian_determinant))
            / 6.0L;
        integrated_power += volume * nonnegative;
        result.electric_field.push_back(field);
        result.current_density.push_back(flux.flux);
        result.joule_density.push_back(density);
    }
    result.total_joule_power_watts = checked_double(
        integrated_power, "Electrothermal 3D total Joule power overflows.");

    result.heated_thermal_problem = problem.thermal;
    auto& sources = result.heated_thermal_problem.cell_source;
    if (!sources.empty() && sources.size() != fluxes.size()) {
        throw std::invalid_argument(
            "Electrothermal 3D thermal cell source count must equal cell count.");
    }
    if (sources.empty()) sources.assign(fluxes.size(), 0.0);
    for (std::size_t i = 0; i < sources.size(); ++i) {
        if (!std::isfinite(sources[i])) {
            throw std::invalid_argument(
                "Electrothermal 3D base thermal cell sources must be finite.");
        }
        sources[i] = checked_double(
            static_cast<long double>(sources[i]) +
            static_cast<long double>(result.joule_density[i]),
            "Electrothermal 3D combined thermal source overflows.");
    }

    result.thermal_solution = fem::solve_scalar_diffusion_3d(
        result.heated_thermal_problem, thermal_options);
    if (!result.thermal_solution.linear_solver.converged) {
        throw std::runtime_error("Electrothermal 3D thermal solve did not converge.");
    }
    result.thermal_balance = fem::evaluate_diffusion_energy_balance_3d(
        result.heated_thermal_problem, result.thermal_solution.values);
    return result;
}

} // namespace nova::physics
