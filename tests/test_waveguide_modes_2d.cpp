#include "test_framework.hpp"
#include "nova/fem/mesh2d.hpp"
#include "nova/physics/waveguide_modes_2d.hpp"
#include "nova/linalg/sparse_symmetric_eigen.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace {

using nova::physics::WaveguideMaterial2D;
using nova::physics::WaveguideModeProblem2D;
using nova::physics::WaveguidePolarization2D;

constexpr double wavelength = 1.55e-6;
constexpr double core_index = 1.50;
constexpr double cladding_index = 1.45;
constexpr double core_half_width = 1.20e-6;
constexpr double domain_half_height = 6.0e-6;

[[nodiscard]] WaveguideModeProblem2D make_slab_problem(
    const std::size_t elements_y,
    const WaveguidePolarization2D polarization) {
    auto rectangle = nova::fem::make_rectangle_mesh_2d(
        -2.0e-6, 2.0e-6, -domain_half_height, domain_half_height,
        2U, elements_y);

    WaveguideModeProblem2D problem;
    problem.mesh = std::move(rectangle.mesh);
    problem.vacuum_wavelength_m = wavelength;
    problem.cladding_index = cladding_index;
    problem.polarization = polarization;

    for (const std::size_t node : rectangle.boundary_nodes) {
        const double y = problem.mesh.nodes[node].y;
        if (std::abs(std::abs(y) - domain_half_height) < 1e-14) {
            problem.dirichlet_boundary_nodes.push_back(node);
        }
    }

    problem.materials.reserve(problem.mesh.elements.size());
    for (const auto& element : problem.mesh.elements) {
        const double center_y = (
            problem.mesh.nodes[element[0]].y +
            problem.mesh.nodes[element[1]].y +
            problem.mesh.nodes[element[2]].y) / 3.0;
        const double index = std::abs(center_y) < core_half_width
            ? core_index : cladding_index;
        problem.materials.push_back({index * index, 1.0});
    }
    return problem;
}

[[nodiscard]] double analytic_slab_fundamental_index(
    const WaveguidePolarization2D polarization,
    const double core_permittivity = core_index * core_index,
    const double core_permeability = 1.0,
    const double cladding_permittivity = cladding_index * cladding_index,
    const double cladding_permeability = 1.0) {
    constexpr double pi = std::numbers::pi_v<double>;
    const double k0 = 2.0 * pi / wavelength;
    const double half_thickness = core_half_width;
    const double core_material_index =
        std::sqrt(core_permittivity * core_permeability);
    const double cladding_material_index =
        std::sqrt(cladding_permittivity * cladding_permeability);
    const double v = k0 * half_thickness *
        std::sqrt(core_material_index * core_material_index -
                  cladding_material_index * cladding_material_index);
    const double q = polarization == WaveguidePolarization2D::TE_Ez
        ? core_permeability / cladding_permeability
        : core_permittivity / cladding_permittivity;

    // Independent bisection of u - atan(q*w/u) = 0 with
    // u=V*cos(theta), w=V*sin(theta); this is the even fundamental mode.
    double lower = 0.0;
    double upper = pi / 2.0;
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double middle = lower + 0.5 * (upper - lower);
        const double residual = v * std::cos(middle) -
            std::atan2(q * std::sin(middle), std::cos(middle));
        if (residual > 0.0) lower = middle;
        else upper = middle;
    }
    const double theta = lower + 0.5 * (upper - lower);
    const double w = v * std::sin(theta);
    const double decay_over_k0 = w / (k0 * half_thickness);
    return std::sqrt(cladding_material_index * cladding_material_index +
                     decay_over_k0 * decay_over_k0);
}

} // namespace

NOVA_TEST(waveguide_te_tm_fundamental_modes_match_independent_slab_reference) {
    for (const auto polarization : {
             WaveguidePolarization2D::TE_Ez,
             WaveguidePolarization2D::TM_Hz}) {
        const auto coarse_problem = make_slab_problem(20U, polarization);
        const auto fine_problem = make_slab_problem(40U, polarization);
        nova::linalg::SparseSymmetricEigenOptions eigen_options;
        eigen_options.max_subspace = 128U;
        const auto coarse = nova::physics::solve_waveguide_modes_2d(
            coarse_problem, 4U, eigen_options);
        const auto fine = nova::physics::solve_waveguide_modes_2d(
            fine_problem, 4U, eigen_options);

        NOVA_REQUIRE(coarse.converged);
        NOVA_REQUIRE(fine.converged);
        NOVA_REQUIRE(!coarse.modes.empty());
        NOVA_REQUIRE(!fine.modes.empty());

        const double reference = analytic_slab_fundamental_index(polarization);
        const double coarse_error = std::abs(coarse.modes.front().effective_index - reference);
        const double fine_error = std::abs(fine.modes.front().effective_index - reference);
        NOVA_REQUIRE(fine_error < coarse_error);
        NOVA_REQUIRE(fine_error < 2.0e-4);

        const auto& fundamental = fine.modes.front();
        NOVA_REQUIRE(fundamental.scalar_field.size() == fine_problem.mesh.nodes.size());
        NOVA_REQUIRE(fundamental.relative_residual < 1.0e-10);
        NOVA_REQUIRE(fundamental.core_confinement_factor > 0.8);
        NOVA_REQUIRE(fundamental.core_confinement_factor <= 1.0);
        NOVA_REQUIRE(fundamental.effective_index > cladding_index);
        NOVA_REQUIRE(fundamental.effective_index < core_index);
        NOVA_REQUIRE_NEAR(
            fundamental.propagation_constant_per_m *
                fundamental.propagation_constant_per_m /
                fundamental.propagation_constant_squared_per_m2,
            1.0,
            1.0e-12);

        for (const std::size_t node : fine_problem.dirichlet_boundary_nodes) {
            NOVA_REQUIRE_NEAR(fundamental.scalar_field[node], 0.0, 0.0);
        }
        // The x-min edge is a natural Neumann boundary, not a constrained edge.
        bool found_unconstrained_side_value = false;
        for (std::size_t node = 0; node < fine_problem.mesh.nodes.size(); ++node) {
            const auto& point = fine_problem.mesh.nodes[node];
            if (std::abs(point.x + 2.0e-6) < 1e-14 &&
                std::abs(point.y) < 1e-14) {
                NOVA_REQUIRE(std::abs(fundamental.scalar_field[node]) > 1e-12);
                found_unconstrained_side_value = true;
            }
        }
        NOVA_REQUIRE(found_unconstrained_side_value);

        for (std::size_t mode = 1; mode < fine.modes.size(); ++mode) {
            NOVA_REQUIRE(
                fine.modes[mode - 1U].effective_index >=
                fine.modes[mode].effective_index);
            NOVA_REQUIRE(fine.modes[mode].core_confinement_factor >= 0.0);
            NOVA_REQUIRE(fine.modes[mode].core_confinement_factor <= 1.0);
        }

        std::cout << (polarization == WaveguidePolarization2D::TE_Ez ? "TE_Ez" : "TM_Hz")
                  << " slab n_eff=" << fundamental.effective_index
                  << " analytic=" << reference
                  << " error=" << fine_error
                  << " residual=" << fundamental.relative_residual << '\n';
    }
}

NOVA_TEST(waveguide_modes_use_permeability_weighted_te_tm_operators) {
    constexpr double core_mu = 1.2;
    constexpr double cladding_mu = 1.0;
    constexpr double core_epsilon = (core_index * core_index) / core_mu;
    constexpr double cladding_epsilon = cladding_index * cladding_index;
    auto problem = make_slab_problem(40U, WaveguidePolarization2D::TE_Ez);
    for (std::size_t index = 0; index < problem.mesh.elements.size(); ++index) {
        const auto& element = problem.mesh.elements[index];
        const double center_y = (
            problem.mesh.nodes[element[0]].y +
            problem.mesh.nodes[element[1]].y +
            problem.mesh.nodes[element[2]].y) / 3.0;
        if (std::abs(center_y) < core_half_width) {
            problem.materials[index] = {core_epsilon, core_mu};
        } else {
            problem.materials[index] = {cladding_epsilon, cladding_mu};
        }
    }

    nova::linalg::SparseSymmetricEigenOptions eigen_options;
    eigen_options.max_subspace = 128U;
    for (const auto polarization : {
             WaveguidePolarization2D::TE_Ez,
             WaveguidePolarization2D::TM_Hz}) {
        problem.polarization = polarization;
        const auto result = nova::physics::solve_waveguide_modes_2d(
            problem, 4U, eigen_options);
        const double reference = analytic_slab_fundamental_index(
            polarization, core_epsilon, core_mu,
            cladding_epsilon, cladding_mu);
        NOVA_REQUIRE(result.converged);
        NOVA_REQUIRE(!result.modes.empty());
        NOVA_REQUIRE_NEAR(result.modes.front().effective_index, reference, 2.0e-4);
        NOVA_REQUIRE(result.modes.front().core_confinement_factor > 0.8);
    }
}

NOVA_TEST(waveguide_modes_are_invariant_to_triangle_orientation) {
    auto problem = make_slab_problem(20U, WaveguidePolarization2D::TE_Ez);
    const auto baseline = nova::physics::solve_waveguide_modes_2d(problem, 2U);
    for (auto& element : problem.mesh.elements) {
        std::swap(element[1], element[2]);
    }
    const auto reversed = nova::physics::solve_waveguide_modes_2d(problem, 2U);
    NOVA_REQUIRE(baseline.converged && reversed.converged);
    NOVA_REQUIRE(!baseline.modes.empty() && !reversed.modes.empty());
    NOVA_REQUIRE_NEAR(
        baseline.modes.front().effective_index,
        reversed.modes.front().effective_index,
        1.0e-11);
}

NOVA_TEST(waveguide_uniform_cladding_has_no_bound_modes) {
    auto problem = make_slab_problem(16U, WaveguidePolarization2D::TE_Ez);
    for (auto& material : problem.materials) {
        material = {cladding_index * cladding_index, 1.0};
    }
    const auto result = nova::physics::solve_waveguide_modes_2d(problem, 3U);
    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.modes.empty());
}

NOVA_TEST(waveguide_modes_reject_invalid_material_mesh_and_constraints) {
    auto problem = make_slab_problem(12U, WaveguidePolarization2D::TE_Ez);

    auto invalid = problem;
    invalid.materials.pop_back();
    bool rejected = false;
    try { (void)nova::physics::solve_waveguide_modes_2d(invalid); }
    catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    invalid = problem;
    invalid.materials.front().relative_permittivity = -1.0;
    rejected = false;
    try { (void)nova::physics::solve_waveguide_modes_2d(invalid); }
    catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    invalid = problem;
    invalid.dirichlet_boundary_nodes.push_back(
        (invalid.mesh.nodes.size() / 2U));
    rejected = false;
    try { (void)nova::physics::solve_waveguide_modes_2d(invalid); }
    catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    invalid = problem;
    invalid.dirichlet_boundary_nodes.push_back(invalid.dirichlet_boundary_nodes.front());
    rejected = false;
    try { (void)nova::physics::solve_waveguide_modes_2d(invalid); }
    catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    invalid = problem;
    invalid.cladding_index = 1.6;
    rejected = false;
    try { (void)nova::physics::solve_waveguide_modes_2d(invalid); }
    catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    rejected = false;
    try { (void)nova::physics::solve_waveguide_modes_2d(problem, 0U); }
    catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    auto invalid_options = nova::linalg::SparseSymmetricEigenOptions{};
    invalid_options.relative_tolerance = 0.0;
    rejected = false;
    try { (void)nova::physics::solve_waveguide_modes_2d(problem, 2U, invalid_options); }
    catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);
}
