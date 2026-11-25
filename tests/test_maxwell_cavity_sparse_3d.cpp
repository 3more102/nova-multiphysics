#include "test_framework.hpp"

#include "nova/physics/maxwell_cavity_sparse_3d.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

using namespace nova::physics;
using namespace nova::linalg;

NOVA_TEST(maxwell_sparse_shift_invert_matches_dense_pec_cube_modes) {
    for (const std::size_t divisions : {2U, 3U}) {
        const auto mesh = nova::fem::make_box_mesh_3d(
            0, 1, 0, 1, 0, 1, divisions, divisions, divisions).mesh;
        const std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size());
        const auto matrices = assemble_maxwell_cavity_3d(mesh, materials);
        const auto dense = solve_maxwell_cavity_3d(matrices, 3);

        MaxwellCavitySparseOptions3D options;
        options.max_modes = 3;
        options.max_subspace = 64;
        options.relative_tolerance = 2.0e-7;
        options.shift_solve_relative_tolerance = 1.0e-12;
        options.projection_relative_tolerance = 1.0e-12;
        const auto sparse = solve_maxwell_cavity_shift_invert_3d(
            mesh, matrices, options);

        NOVA_REQUIRE(dense.converged);
        NOVA_REQUIRE(sparse.converged);
        NOVA_REQUIRE(sparse.modes.size() == 3U);
        NOVA_REQUIRE(sparse.gradient_nullspace_dimension ==
                     (divisions - 1U) * (divisions - 1U) * (divisions - 1U));
        NOVA_REQUIRE(sparse.operator_applications == sparse.subspace_dimension);
        NOVA_REQUIRE(sparse.subspace_dimension <= options.max_subspace);
        NOVA_REQUIRE(sparse.shift > 0.0);
        NOVA_REQUIRE(sparse.shift_solve_iterations > 0U);

        std::vector<std::vector<double>> reduced_modes;
        for (std::size_t mode_index = 0; mode_index < sparse.modes.size(); ++mode_index) {
            const auto& mode = sparse.modes[mode_index];
            NOVA_REQUIRE_NEAR(
                mode.vacuum_wavenumber_squared,
                dense.modes[mode_index].vacuum_wavenumber_squared,
                3.0e-6);
            NOVA_REQUIRE(mode.relative_residual < options.relative_tolerance);
            NOVA_REQUIRE(mode.frequency_hz > 0.0);
            NOVA_REQUIRE(mode.edge_field.size() == matrices.edges.size());

            std::vector<double> reduced(matrices.free_edges.size(), 0.0);
            for (std::size_t i = 0; i < matrices.free_edges.size(); ++i) {
                reduced[i] = mode.edge_field[matrices.free_edges[i]];
            }
            const auto mass_vector = matrices.mass.multiply(reduced);
            double mass_norm = 0.0;
            for (std::size_t i = 0; i < reduced.size(); ++i) {
                mass_norm += reduced[i] * mass_vector[i];
            }
            NOVA_REQUIRE_NEAR(mass_norm, 1.0, 2.0e-9);
            for (const auto& previous : reduced_modes) {
                const auto previous_mass = matrices.mass.multiply(previous);
                double inner = 0.0;
                for (std::size_t i = 0; i < reduced.size(); ++i) {
                    inner += reduced[i] * previous_mass[i];
                }
                NOVA_REQUIRE_NEAR(inner, 0.0, 2.0e-8);
            }

            std::vector<unsigned char> boundary(mesh.nodes.size(), 0U);
            for (const auto& face : nova::fem::extract_boundary_faces_3d(mesh)) {
                for (const auto node : face) boundary[node] = 1U;
            }
            std::vector<double> gradient_residual(mesh.nodes.size(), 0.0);
            for (std::size_t i = 0; i < matrices.free_edges.size(); ++i) {
                const auto& edge = matrices.edges[matrices.free_edges[i]];
                if (boundary[edge[0]] == 0U) gradient_residual[edge[0]] -= mass_vector[i];
                if (boundary[edge[1]] == 0U) gradient_residual[edge[1]] += mass_vector[i];
            }
            double gradient_norm = 0.0;
            for (const double value : gradient_residual) {
                gradient_norm = std::hypot(gradient_norm, value);
            }
            NOVA_REQUIRE(gradient_norm < 2.0e-8);

            for (std::size_t edge = 0; edge < matrices.edges.size(); ++edge) {
                if (!std::binary_search(
                        matrices.free_edges.begin(), matrices.free_edges.end(), edge)) {
                    NOVA_REQUIRE(mode.edge_field[edge] == 0.0);
                }
            }
            reduced_modes.push_back(std::move(reduced));
        }
    }
}

NOVA_TEST(maxwell_sparse_shift_invert_preserves_dielectric_frequency_scaling) {
    const auto mesh = nova::fem::make_box_mesh_3d(
        0, 1, 0, 1, 0, 1, 2, 2, 2).mesh;
    std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size());

    MaxwellCavitySparseOptions3D options;
    options.max_modes = 1;
    options.max_subspace = 32;
    const auto vacuum = solve_maxwell_cavity_shift_invert_3d(
        mesh, assemble_maxwell_cavity_3d(mesh, materials), options);

    for (auto& material : materials) material.relative_permittivity = 4.0;
    const auto dielectric = solve_maxwell_cavity_shift_invert_3d(
        mesh, assemble_maxwell_cavity_3d(mesh, materials), options);

    NOVA_REQUIRE(vacuum.converged && dielectric.converged);
    NOVA_REQUIRE_NEAR(
        dielectric.modes[0].vacuum_wavenumber_squared /
            vacuum.modes[0].vacuum_wavenumber_squared,
        0.25, 2.0e-8);
    NOVA_REQUIRE_NEAR(
        dielectric.modes[0].frequency_hz / vacuum.modes[0].frequency_hz,
        0.5, 2.0e-8);
}

NOVA_TEST(maxwell_sparse_shift_invert_rejects_wrong_mesh_and_invalid_options) {
    const auto mesh_a = nova::fem::make_box_mesh_3d(
        0, 1, 0, 1, 0, 1, 2, 2, 2).mesh;
    const auto mesh_b = nova::fem::make_box_mesh_3d(
        0, 1, 0, 1, 0, 1, 3, 3, 3).mesh;
    const std::vector<MaxwellCellMaterial3D> materials_a(mesh_a.elements.size());
    const auto matrices_a = assemble_maxwell_cavity_3d(mesh_a, materials_a);

    bool rejected = false;
    try {
        (void)solve_maxwell_cavity_shift_invert_3d(mesh_b, matrices_a);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);

    MaxwellCavitySparseOptions3D invalid;
    invalid.max_subspace = 513;
    rejected = false;
    try {
        (void)solve_maxwell_cavity_shift_invert_3d(mesh_a, matrices_a, invalid);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}
