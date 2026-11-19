#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/poisson_1d.hpp"
#include "nova/linalg/csr_matrix.hpp"
#include "nova/physics/electrostatics_3d.hpp"
#include "nova/io/gmsh_tetrahedral.hpp"

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
    return 0;
}
