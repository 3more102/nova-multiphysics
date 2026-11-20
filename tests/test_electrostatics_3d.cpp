#include "test_framework.hpp"
#include "nova/physics/electrostatics_3d.hpp"
#include "nova/physics/electrode_charge_3d.hpp"

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
    p.free_charge_density = [](double, double, double) { return 2.0; };
    p.cell_free_charge_density = {2.0};
    p.potential_constraints = {{0, 0.0}, {2, 0.0}, {3, 0.0}};

    // A single free P1 DOF: K_11 = epsilon*volume = 1/6;
    // F_1 = rho*volume/4 = 1/6, so V_1 = F_1/K_11 = 1 V.
    const auto result = solve_electrostatics_3d(p);
    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE_NEAR(result.potential[1], 1.0, 1e-11);
    NOVA_REQUIRE_NEAR(result.electric_field[0].x, -1.0, 1e-11);
    NOVA_REQUIRE_NEAR(result.stored_energy_joules, 1.0 / 12.0, 1e-11);
}


NOVA_TEST(electrostatics_3d_cell_only_charge_source) {
    Electrostatics3DProblem p;
    p.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    p.cell_permittivity = {1.0};
    p.cell_free_charge_density = {4.0};
    p.potential_constraints = {{0, 0.0}, {2, 0.0}, {3, 0.0}};

    // Missing spatial charge means zero; the cell source alone gives V_1 = 1 V.
    const auto result = solve_electrostatics_3d(p);
    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE_NEAR(result.potential[1], 1.0, 1e-11);
    NOVA_REQUIRE_NEAR(result.stored_energy_joules, 1.0 / 12.0, 1e-11);
}

NOVA_TEST(electrostatics_3d_prescribed_displacement_recovers_pure_neumann_affine_field) {
    using namespace nova::fem;
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 2U, 2U, 2U);
    Electrostatics3DProblem p;
    p.mesh = box.mesh;
    p.permittivity = [](double, double, double) { return 2.0; };
    for (const auto& face : box.boundary_faces) {
        const auto on_x = [&](double x) {
            return std::all_of(face.begin(), face.end(),
                [&](std::size_t vertex) { return p.mesh.nodes[vertex].x == x; });
        };
        const double outward_displacement = on_x(0.0) ? 2.0 : (on_x(1.0) ? -2.0 : 0.0);
        p.displacement_flux_faces.push_back({face,
            [outward_displacement](double, double, double) { return outward_displacement; }});
    }

    // Pure-Neumann compatibility fixes V only up to a constant; NOVA returns
    // the zero-volume-mean representative V=x-1/2.
    const auto result = solve_electrostatics_3d(p);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t i = 0; i < result.potential.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.potential[i], p.mesh.nodes[i].x - 0.5, 1e-9);
    }
    for (const auto& field : result.electric_field) {
        NOVA_REQUIRE_NEAR(field.x, -1.0, 1e-9);
        NOVA_REQUIRE_NEAR(field.y, 0.0, 1e-9);
        NOVA_REQUIRE_NEAR(field.z, 0.0, 1e-9);
    }
    for (const auto& displacement : result.displacement) {
        NOVA_REQUIRE_NEAR(displacement.x, -2.0, 1e-9);
        NOVA_REQUIRE_NEAR(displacement.y, 0.0, 1e-9);
        NOVA_REQUIRE_NEAR(displacement.z, 0.0, 1e-9);
    }
    NOVA_REQUIRE_NEAR(result.stored_energy_joules, 1.0, 1e-9);
}

NOVA_TEST(electrostatics_3d_electrode_charge_parallel_plate) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2, 2, 2
    );
    Electrostatics3DProblem p;
    p.mesh = box.mesh;
    p.permittivity = [](double, double, double) { return 2.0; };
    std::vector<nova::fem::BoundaryFace3D> electrodes;
    std::vector<int> ids;
    for (const auto& face : box.boundary_faces) {
        if (on_x(p.mesh, face, 0.0) || on_x(p.mesh, face, 1.0)) {
            const bool high = on_x(p.mesh, face, 1.0);
            p.potential_faces.push_back({
                face, [high](double, double, double) { return high ? 10.0 : 0.0; }
            });
            // Deliberately reverse the input winding: the physical sign
            // must follow the dielectric's outward normal, not this order.
            electrodes.push_back({face[2], face[1], face[0]});
            ids.push_back(high ? 20 : 10);
        }
    }

    const auto solved = solve_electrostatics_3d(p);
    const auto charges = nova::physics::evaluate_electrode_charges_3d(
        p, solved, electrodes, ids
    );
    NOVA_REQUIRE(charges.size() == 2U);
    NOVA_REQUIRE(charges[0].electrode_id == 10);
    NOVA_REQUIRE(charges[1].electrode_id == 20);
    NOVA_REQUIRE_NEAR(charges[0].area_square_metres, 1.0, 1e-12);
    NOVA_REQUIRE_NEAR(charges[1].area_square_metres, 1.0, 1e-12);
    NOVA_REQUIRE_NEAR(charges[0].potential_volts, 0.0, 1e-12);
    NOVA_REQUIRE_NEAR(charges[1].potential_volts, 10.0, 1e-12);
    NOVA_REQUIRE_NEAR(charges[0].conductor_charge_coulombs, -20.0, 1e-8);
    NOVA_REQUIRE_NEAR(charges[1].conductor_charge_coulombs, 20.0, 1e-8);
    NOVA_REQUIRE_NEAR(charges[0].dielectric_outward_charge_coulombs, 20.0, 1e-8);
    NOVA_REQUIRE_NEAR(charges[1].dielectric_outward_charge_coulombs, -20.0, 1e-8);
    // Here there is no volume charge or other electrode. C=Q_high/dV=2 F;
    // the independently integrated field energy must equal 0.5*C*dV^2.
    const double capacitance = charges[1].conductor_charge_coulombs /
        (charges[1].potential_volts - charges[0].potential_volts);
    NOVA_REQUIRE_NEAR(capacitance, 2.0, 1e-9);
    NOVA_REQUIRE_NEAR(solved.stored_energy_joules, 100.0, 1e-7);
}

NOVA_TEST(electrostatics_3d_electrode_charge_rejects_invalid_face_sets) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1
    );
    Electrostatics3DProblem p;
    p.mesh = box.mesh;
    Electrostatics3DResult solution;
    for (std::size_t i = 0; i < p.mesh.nodes.size(); ++i) {
        p.potential_constraints.push_back({i, p.mesh.nodes[i].x});
        solution.potential.push_back(p.mesh.nodes[i].x);
    }
    solution.displacement.assign(p.mesh.elements.size(), {-1.0, 0.0, 0.0});

    nova::fem::BoundaryFace3D low{}, high{};
    bool found_low = false, found_high = false;
    for (const auto& face : box.boundary_faces) {
        if (on_x(p.mesh, face, 0.0) && !found_low) {
            low = face;
            found_low = true;
        }
        if (on_x(p.mesh, face, 1.0) && !found_high) {
            high = face;
            found_high = true;
        }
    }
    NOVA_REQUIRE(found_low && found_high);
    const auto reject = [](const auto& operation) {
        bool invalid = false;
        try { operation(); }
        catch (const std::invalid_argument&) { invalid = true; }
        NOVA_REQUIRE(invalid);
    };

    // One group cannot combine two different fixed potentials.
    reject([&] {
        (void)nova::physics::evaluate_electrode_charges_3d(
            p, solution, std::vector<nova::fem::BoundaryFace3D>{low, high},
            std::vector<int>{7, 7}
        );
    });
    // Faces must be unique, even if their tags differ.
    reject([&] {
        (void)nova::physics::evaluate_electrode_charges_3d(
            p, solution, std::vector<nova::fem::BoundaryFace3D>{low, low},
            std::vector<int>{7, 8}
        );
    });
    // A selected face with unconstrained vertices is not a conductor electrode.
    p.potential_constraints.clear();
    p.potential_constraints.push_back({low[0], 0.0});
    reject([&] {
        (void)nova::physics::evaluate_electrode_charges_3d(
            p, solution, std::vector<nova::fem::BoundaryFace3D>{low},
            std::vector<int>{7}
        );
    });
    // NaN/Inf displacement is never silently reported as a finite charge.
    p.potential_constraints.clear();
    for (std::size_t i = 0; i < p.mesh.nodes.size(); ++i) {
        p.potential_constraints.push_back({i, p.mesh.nodes[i].x});
    }
    solution.displacement[0].x = std::numeric_limits<double>::infinity();
    reject([&] {
        (void)nova::physics::evaluate_electrode_charges_3d(
            p, solution, std::vector<nova::fem::BoundaryFace3D>{low},
            std::vector<int>{7}
        );
    });
}
