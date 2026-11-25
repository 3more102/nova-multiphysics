#include "nova/physics/maxwell_cavity_sparse_3d.hpp"

#include <iomanip>
#include <iostream>
#include <vector>

int main() {
    using namespace nova::physics;
    const auto mesh = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 3, 3, 3).mesh;
    const std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size());
    const auto matrices = assemble_maxwell_cavity_3d(mesh, materials);

    MaxwellCavitySparseOptions3D options;
    options.max_modes = 3;
    options.max_subspace = 64;
    const auto result = solve_maxwell_cavity_shift_invert_3d(
        mesh, matrices, options);

    std::cout << std::setprecision(12)
              << "NOVA sparse PEC cavity shift-invert\n"
              << "Tetrahedra: " << mesh.elements.size() << '\n'
              << "Free edge DOFs: " << matrices.free_edges.size() << '\n'
              << "Exact gradient nullity: " << result.gradient_nullspace_dimension << '\n'
              << "Observed additional zero modes: " << result.additional_zero_modes_detected << '\n'
              << "Krylov subspace: " << result.subspace_dimension << '\n'
              << "Positive shift: " << result.shift << '\n'
              << "Shift-solve CG iterations: " << result.shift_solve_iterations << '\n'
              << "Projection CG iterations: " << result.projection_solve_iterations << '\n'
              << "Converged: " << (result.converged ? "yes" : "no") << '\n';

    for (std::size_t i = 0; i < result.modes.size(); ++i) {
        const auto& mode = result.modes[i];
        std::cout << "Mode " << (i + 1U)
                  << ": lambda=" << mode.vacuum_wavenumber_squared
                  << " frequency_Hz=" << mode.frequency_hz
                  << " relative_residual=" << mode.relative_residual << '\n';
    }
    return result.converged ? 0 : 1;
}
