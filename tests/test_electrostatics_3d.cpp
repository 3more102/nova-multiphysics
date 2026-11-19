#include "test_framework.hpp"
#include "nova/physics/electrostatics_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

using nova::physics::Electrostatics3DProblem;
using nova::physics::solve_electrostatics_3d;

bool on_x(const nova::fem::TetrahedralMesh3D& mesh,
          const nova::fem::BoundaryFace3D& face, double x) {
    return std::all_of(face.begin(), face.end(),
        [&](std::size_t node) { return mesh.nodes[node].x == x; });
}

} // namespace

NOVA_TEST(electrostatics_3d_parallel_plate_stored_energy) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2, 2, 2
    );
    Electrostatics3DProblem p;
    p.mesh = box.mesh;
    p.permittivity = [](double, double, double) { return 2.0; };
    for (const auto& face : box.boundary_faces) {
        if (on_x(p.mesh, face, 0.0)) {
            p.potential_faces.push_back({
                face, [](double, double, double) { return 0.0; }
            });
        } else if (on_x(p.mesh, face, 1.0)) {
            p.potential_faces.push_back({
                face, [](double, double, double) { return 10.0; }
            });
        }
    }

    const auto result = solve_electrostatics_3d(p);
    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE(result.potential.size() == p.mesh.nodes.size());
    NOVA_REQUIRE(result.electric_field.size() == p.mesh.elements.size());
    for (std::size_t i = 0; i < p.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.potential[i], 10.0 * p.mesh.nodes[i].x, 1e-8);
    }
    for (std::size_t i = 0; i < result.electric_field.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.electric_field[i].x, -10.0, 1e-8);
        NOVA_REQUIRE_NEAR(result.electric_field[i].y, 0.0, 1e-8);
        NOVA_REQUIRE_NEAR(result.electric_field[i].z, 0.0, 1e-8);
        NOVA_REQUIRE_NEAR(result.displacement[i].x, -20.0, 1e-8);
    }
    // C = epsilon A / gap = 2 F; U = 1/2 C (10 V)^2 = 100 J.
    NOVA_REQUIRE_NEAR(result.stored_energy_joules, 100.0, 1e-7);
}

NOVA_TEST(electrostatics_3d_anisotropic_tensor_affine_field) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1
    );
    Electrostatics3DProblem p;
    p.mesh = box.mesh;
    p.cell_permittivity_tensors.assign(
        p.mesh.elements.size(), {4.0, 5.0, 6.0, 0.25, -0.5, 0.75}
    );
    for (const auto& face : box.boundary_faces) {
        p.potential_faces.push_back({
            face,
            [](double x, double y, double z) {
                return x + 2.0 * y - 3.0 * z;
            }
        });
    }

    const auto result = solve_electrostatics_3d(p);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t i = 0; i < p.mesh.nodes.size(); ++i) {
        const auto& node = p.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(result.potential[i], node.x + 2.0 * node.y - 3.0 * node.z, 1e-10);
    }
    for (std::size_t i = 0; i < result.electric_field.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.electric_field[i].x, -1.0, 1e-10);
        NOVA_REQUIRE_NEAR(result.electric_field[i].y, -2.0, 1e-10);
        NOVA_REQUIRE_NEAR(result.electric_field[i].z, 3.0, 1e-10);
        NOVA_REQUIRE_NEAR(result.displacement[i].x, -6.0, 1e-10);
        NOVA_REQUIRE_NEAR(result.displacement[i].y, -8.0, 1e-10);
        NOVA_REQUIRE_NEAR(result.displacement[i].z, 17.0, 1e-10);
    }
    // E dot D = 73 J/m^3 in this SI-material benchmark.
    NOVA_REQUIRE_NEAR(result.stored_energy_joules, 36.5, 1e-10);
}

NOVA_TEST(electrostatics_3d_rejects_nonphysical_inputs) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1
    );
    Electrostatics3DProblem p;
    p.mesh = box.mesh;
    p.permittivity = [](double, double, double) { return -1.0; };
    p.potential_constraints.push_back({0, 0.0});
    bool negative_rejected = false;
    try { (void)solve_electrostatics_3d(p); }
    catch (const std::invalid_argument&) { negative_rejected = true; }
    NOVA_REQUIRE(negative_rejected);

    p.permittivity = [](double, double, double) { return 1.0; };
    p.cell_permittivity.assign(p.mesh.elements.size(), 2.0);
    bool double_material_rejected = false;
    try { (void)solve_electrostatics_3d(p); }
    catch (const std::invalid_argument&) { double_material_rejected = true; }
    NOVA_REQUIRE(double_material_rejected);

    p.cell_permittivity.clear();
    p.free_charge_density = [](double, double, double) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    bool charge_rejected = false;
    try { (void)solve_electrostatics_3d(p); }
    catch (const std::invalid_argument&) { charge_rejected = true; }
    NOVA_REQUIRE(charge_rejected);
}

NOVA_TEST(electrostatics_3d_unit_tetra_space_charge_exact_fem) {
    Electrostatics3DProblem p;
    p.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    p.cell_permittivity = {1.0};
    p.free_charge_density = [](double, double, double) { return 4.0; };
    p.potential_constraints = {{0, 0.0}, {2, 0.0}, {3, 0.0}};

    // A single free P1 DOF: K_11 = epsilon*volume = 1/6;
    // F_1 = rho*volume/4 = 1/6, so V_1 = F_1/K_11 = 1 V.
    const auto result = solve_electrostatics_3d(p);
    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE_NEAR(result.potential[1], 1.0, 1e-11);
    NOVA_REQUIRE_NEAR(result.electric_field[0].x, -1.0, 1e-11);
    NOVA_REQUIRE_NEAR(result.stored_energy_joules, 1.0 / 12.0, 1e-11);
}
