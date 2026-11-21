#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/poisson_1d.hpp"
#include "nova/linalg/csr_matrix.hpp"
#include "nova/linalg/conjugate_gradient.hpp"
#include "nova/linalg/incomplete_cholesky.hpp"
#include "nova/physics/electrostatics_3d.hpp"
#include "nova/physics/electrode_charge_3d.hpp"
#include "nova/io/gmsh_tetrahedral.hpp"
#include "nova/io/thermal_checkpoint_3d.hpp"

#include <array>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

int main() {
    // Exercise 1D FEM and CSR algebra across the installed library boundary.
    nova::fem::Poisson1DProblem p;
    p.elements = 16;
    p.source = [](double) { return 2.0; };
    const auto one_d = nova::fem::solve_poisson_1d(p);
    if (one_d.values.size() != 17 ||
        std::abs(one_d.values[8] - 0.25) > 1.0e-7) return 1;

    const auto a = nova::linalg::CsrMatrix::from_triplets(2, 2,
        {{0, 0, 2.0}, {1, 1, 3.0}});
    const std::vector<double> x{4.0, 5.0};
    const auto y = a.multiply(x);
    if (y.size() != 2 || std::abs(y[0] - 8.0) > 1.0e-12 ||
        std::abs(y[1] - 15.0) > 1.0e-12) return 2;

    // Exercise the installed IC(0) factor and CG selection API.
    const auto spd = nova::linalg::CsrMatrix::from_triplets(2, 2,
        {{0, 0, 4.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 3.0}});
    const std::array<double, 2> rhs{1.0, 2.0};
    nova::linalg::SolverOptions ic_options;
    ic_options.preconditioner = nova::linalg::PreconditionerType::IncompleteCholesky0;
    ic_options.relative_tolerance = 1.0e-12;
    ic_options.absolute_tolerance = 1.0e-14;
    const auto ic_solution =
        nova::linalg::solve_conjugate_gradient(spd, rhs, ic_options);
    if (!ic_solution.converged || ic_solution.iterations != 1 ||
        std::abs(ic_solution.solution[0] - 1.0 / 11.0) > 1.0e-12 ||
        std::abs(ic_solution.solution[1] - 7.0 / 11.0) > 1.0e-12) return 8;

    nova::linalg::CachedConjugateGradient cached_cg;
    const std::array<double, 2> second_exact{-1.0, 2.0};
    const auto cached_first = cached_cg.solve(spd, rhs, ic_options);
    const auto cached_second =
        cached_cg.solve(spd, spd.multiply(second_exact), ic_options);
    if (!cached_first.converged || !cached_second.converged ||
        cached_cg.factor_builds() != 1 ||
        std::abs(cached_second.solution[0] + 1.0) > 1.0e-12 ||
        std::abs(cached_second.solution[1] - 2.0) > 1.0e-12) return 9;

    // Independent installed 3D API consumer: -Delta(u)=0, u=x, with
    // fixed temperatures at x=0 and x=1 and insulated remaining faces.
    using namespace nova::fem;
    const auto box = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    for (const auto& face : box.boundary_faces) {
        const bool left = problem.mesh.nodes[face[0]].x == 0.0 &&
            problem.mesh.nodes[face[1]].x == 0.0 &&
            problem.mesh.nodes[face[2]].x == 0.0;
        const bool right = problem.mesh.nodes[face[0]].x == 1.0 &&
            problem.mesh.nodes[face[1]].x == 1.0 &&
            problem.mesh.nodes[face[2]].x == 1.0;
        if (left || right) {
            problem.dirichlet_faces.push_back({face,
                [](double px, double, double) { return px; }});
        }
    }
    const auto three_d = solve_scalar_diffusion_3d(problem);
    if (!three_d.linear_solver.converged ||
        three_d.values.size() != problem.mesh.nodes.size()) return 3;
    for (std::size_t i = 0; i < three_d.values.size(); ++i) {
        if (std::abs(three_d.values[i] - problem.mesh.nodes[i].x) > 1.0e-9)
            return 4;
    }
    // Exercise the installed experimental electrostatics API and exported symbol.
    nova::physics::Electrostatics3DProblem electrostatic;
    electrostatic.mesh = box.mesh;
    electrostatic.permittivity = [](double, double, double) { return 2.0; };
    std::vector<nova::fem::BoundaryFace3D> electrode_faces;
    std::vector<int> electrode_ids;
    for (const auto& face : box.boundary_faces) {
        const bool left = electrostatic.mesh.nodes[face[0]].x == 0.0 &&
            electrostatic.mesh.nodes[face[1]].x == 0.0 &&
            electrostatic.mesh.nodes[face[2]].x == 0.0;
        const bool right = electrostatic.mesh.nodes[face[0]].x == 1.0 &&
            electrostatic.mesh.nodes[face[1]].x == 1.0 &&
            electrostatic.mesh.nodes[face[2]].x == 1.0;
        if (left || right) {
            electrostatic.potential_faces.push_back({face,
                [](double px, double, double) { return 10.0 * px; }});
            electrode_faces.push_back(face);
            electrode_ids.push_back(right ? 2 : 1);
        }
    }
    const auto capacitor = nova::physics::solve_electrostatics_3d(electrostatic);
    if (!capacitor.linear_solver.converged ||
        capacitor.potential.size() != electrostatic.mesh.nodes.size() ||
        std::abs(capacitor.stored_energy_joules - 100.0) > 1.0e-8) return 5;
    for (std::size_t i = 0; i < capacitor.potential.size(); ++i) {
        if (std::abs(capacitor.potential[i] - 10.0 * electrostatic.mesh.nodes[i].x) > 1.0e-9)
            return 6;
    }
    // Verify the installed fixed-electrode charge API and exported symbol.
    const auto electrode_charges = nova::physics::evaluate_electrode_charges_3d(
        electrostatic, capacitor, electrode_faces, electrode_ids);
    if (electrode_charges.size() != 2 ||
        electrode_charges[0].electrode_id != 1 ||
        electrode_charges[1].electrode_id != 2 ||
        std::abs(electrode_charges[0].conductor_charge_coulombs + 20.0) > 1.0e-8 ||
        std::abs(electrode_charges[1].conductor_charge_coulombs - 20.0) > 1.0e-8 ||
        std::abs(electrode_charges[1].area_square_metres - 1.0) > 1.0e-12) return 10;
    const double installed_capacitance =
        electrode_charges[1].conductor_charge_coulombs /
        (electrode_charges[1].potential_volts - electrode_charges[0].potential_volts);
    if (std::abs(installed_capacitance - 2.0) > 1.0e-9) return 11;
    // Verify the newly installed generic Gmsh importer as an SDK consumer.
    // This minimal MSH 4.1 file has sparse-compatible block records and no
    // physical groups; topology and exterior faces are reconstructed by NOVA.
    const std::string msh41 = R"msh($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
3 1 0 4
1 2 3 4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
$EndElements
)msh";
    std::istringstream mesh_stream(msh41);
    const auto imported = nova::io::read_gmsh_tetrahedra(mesh_stream);
    if (imported.mesh.nodes.size() != 4 ||
        imported.mesh.elements.size() != 1 ||
        imported.volume_physical_tags.size() != 1 ||
        imported.volume_physical_tags[0] != 0 ||
        imported.exterior_faces.size() != 4) return 7;

    // Verify the installed multi-terminal capacitance-matrix API and symbol.
    std::vector<nova::physics::ElectrodeBoundaryGroup3D> electrode_groups{
        {1, {}}, {2, {}}
    };
    for (const auto& face : box.boundary_faces) {
        const bool left = electrostatic.mesh.nodes[face[0]].x == 0.0 &&
            electrostatic.mesh.nodes[face[1]].x == 0.0 &&
            electrostatic.mesh.nodes[face[2]].x == 0.0;
        const bool right = electrostatic.mesh.nodes[face[0]].x == 1.0 &&
            electrostatic.mesh.nodes[face[1]].x == 1.0 &&
            electrostatic.mesh.nodes[face[2]].x == 1.0;
        if (left) electrode_groups[0].faces.push_back(face);
        if (right) electrode_groups[1].faces.push_back(face);
    }
    const auto capacitance_matrix =
        nova::physics::extract_capacitance_matrix_3d(electrostatic, electrode_groups);
    if (capacitance_matrix.electrode_ids != std::vector<int>{1, 2} ||
        capacitance_matrix.matrix_farads.size() != 4 ||
        std::abs(capacitance_matrix.matrix_farads[0] - 2.0) > 1.0e-8 ||
        std::abs(capacitance_matrix.matrix_farads[1] + 2.0) > 1.0e-8 ||
        std::abs(capacitance_matrix.matrix_farads[2] + 2.0) > 1.0e-8 ||
        std::abs(capacitance_matrix.matrix_farads[3] - 2.0) > 1.0e-8 ||
        capacitance_matrix.max_reciprocity_defect_farads > 1.0e-8) return 13;

    // Verify the installed SDK exports the versioned thermal checkpoint API.
    const std::vector<double> saved_temperature(box.mesh.nodes.size(), 301.25);
    std::stringstream checkpoint_stream;
    nova::io::write_transient_heat_checkpoint_3d(
        checkpoint_stream, box.mesh, 1.75, saved_temperature);
    std::istringstream checkpoint_input(checkpoint_stream.str());
    const auto restored = nova::io::read_transient_heat_checkpoint_3d(
        checkpoint_input, box.mesh);
    if (restored.time != 1.75 || restored.temperature != saved_temperature) return 12;
    return 0;
}
