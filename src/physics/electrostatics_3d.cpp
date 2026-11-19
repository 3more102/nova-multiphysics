#include "nova/physics/electrostatics_3d.hpp"

#include "nova/fem/p1_flux_3d.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace nova::physics {

Electrostatics3DResult solve_electrostatics_3d(
    const Electrostatics3DProblem& problem,
    const linalg::SolverOptions& options
) {
    fem::ScalarDiffusion3DProblem weak;
    weak.mesh = problem.mesh;
    weak.diffusivity = problem.permittivity;
    weak.diffusivity_tensor = problem.permittivity_tensor;
    weak.cell_diffusivity = problem.cell_permittivity;
    weak.cell_diffusivity_tensors = problem.cell_permittivity_tensors;
    weak.source = problem.free_charge_density
        ? problem.free_charge_density
        : fem::ScalarField3D{[](double, double, double) { return 0.0; }};
    weak.cell_source = problem.cell_free_charge_density;
    weak.dirichlet_constraints = problem.potential_constraints;
    weak.dirichlet_faces = problem.potential_faces;
    weak.neumann_faces = problem.displacement_flux_faces;

    auto solved = fem::solve_scalar_diffusion_3d(weak, options);
    auto fluxes = fem::evaluate_p1_tetrahedral_flux_3d(weak, solved.values);

    Electrostatics3DResult result;
    result.potential = std::move(solved.values);
    result.linear_solver = std::move(solved.linear_solver);
    result.electric_field.reserve(fluxes.size());
    result.displacement.reserve(fluxes.size());

    long double energy = 0.0L;
    for (std::size_t index = 0; index < fluxes.size(); ++index) {
        const auto& flux = fluxes[index];
        const fem::Point3D field{
            -flux.gradient.x, -flux.gradient.y, -flux.gradient.z
        };
        const long double xx = static_cast<long double>(field.x) * flux.flux.x;
        const long double yy = static_cast<long double>(field.y) * flux.flux.y;
        const long double zz = static_cast<long double>(field.z) * flux.flux.z;
        const long double density_twice = xx + yy + zz;

        // Positive-definite epsilon implies E dot D >= 0. Accept only
        // floating-point cancellation at the last few significant bits.
        const long double scale = std::abs(xx) + std::abs(yy) + std::abs(zz);
        if (!std::isfinite(density_twice) ||
            density_twice < -256.0L * std::numeric_limits<long double>::epsilon() * scale) {
            throw std::overflow_error("Electrostatic energy density is non-finite or negative.");
        }
        const auto basis = fem::evaluate_tetrahedron_p1_physical_gradients_3d(
            weak.mesh, index
        );
        const long double volume =
            std::abs(static_cast<long double>(basis.signed_jacobian_determinant)) / 6.0L;
        energy += 0.5L * volume * std::max(0.0L, density_twice);
        result.electric_field.push_back(field);
        result.displacement.push_back(flux.flux);
    }
    if (!std::isfinite(energy) ||
        energy > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Electrostatic stored energy overflows double precision.");
    }
    result.stored_energy_joules = static_cast<double>(energy);
    return result;
}

} // namespace nova::physics
