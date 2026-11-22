#include "test_framework.hpp"

#include "nova/physics/linear_elasticity_3d.hpp"

#include <cmath>
#include <exception>
#include <limits>
#include <set>
#include <utility>
#include <vector>

namespace {

using namespace nova::physics;
using namespace nova::fem;

[[nodiscard]] ElasticVector3D affine_displacement(const Point3D& p) {
    return {0.01*p.x + 0.02*p.y,
            -0.015*p.y + 0.03*p.z,
            0.005*p.z + 0.004*p.x};
}

[[nodiscard]] LinearElasticity3DProblem make_model(
    const StructuredBoxMesh3D& box, double young, double poisson
) {
    LinearElasticity3DProblem problem;
    problem.mesh = box.mesh;
    problem.cell_materials.assign(
        box.mesh.elements.size(), {young, poisson});
    return problem;
}

void fix_exterior(
    LinearElasticity3DProblem& problem,
    const std::vector<BoundaryFace3D>& faces,
    ElasticVector3D (*field)(const Point3D&)
) {
    std::set<std::size_t> nodes;
    for (const auto& face : faces) {
        nodes.insert(face.begin(), face.end());
    }
    for (const auto node : nodes) {
        problem.prescribed_displacements.push_back(
            {node, field(problem.mesh.nodes[node])});
    }
}

template <class Action>
[[nodiscard]] bool rejects(Action&& action) {
    try {
        (void)action();
    } catch (const std::exception&) {
        return true;
    }
    return false;
}

} // namespace

NOVA_TEST(elasticity_3d_affine_constant_strain_patch_and_energy) {
    const auto box = make_box_mesh_3d(
        0.0,1.0,0.0,1.0,0.0,1.0,2U,2U,2U);
    auto problem = make_model(box, 2000.0, 0.25);
    fix_exterior(problem, box.boundary_faces, affine_displacement);

    nova::linalg::SolverOptions options;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    options.relative_tolerance = 1.0e-12;
    const auto solved = solve_linear_elasticity_3d(problem, options);
    NOVA_REQUIRE(solved.linear_solver.converged);
    NOVA_REQUIRE(solved.displacement.size() == problem.mesh.nodes.size());
    NOVA_REQUIRE(solved.cell_strain.size() == problem.mesh.elements.size());
    NOVA_REQUIRE(solved.cell_stress.size() == problem.mesh.elements.size());

    for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
        const auto exact = affine_displacement(problem.mesh.nodes[i]);
        for (std::size_t c = 0; c < 3U; ++c) {
            NOVA_REQUIRE_NEAR(solved.displacement[i][c], exact[c], 1.0e-9);
        }
    }
    const ElasticTensor6 strain{0.01, -0.015, 0.005, 0.01, 0.015, 0.002};
    const ElasticTensor6 stress{16.0, -24.0, 8.0, 16.0, 24.0, 3.2};
    for (std::size_t e = 0; e < solved.cell_strain.size(); ++e) {
        for (std::size_t c = 0; c < 6U; ++c) {
            NOVA_REQUIRE_NEAR(solved.cell_strain[e][c], strain[c], 1.0e-9);
            NOVA_REQUIRE_NEAR(solved.cell_stress[e][c], stress[c], 1.0e-6);
        }
    }
    // E=2000, nu=0.25 => lambda=mu=800; unit volume.
    // W = 0.5 * sigma:epsilon = 0.8064 J.
    NOVA_REQUIRE_NEAR(solved.strain_energy_joules, 0.8064, 1.0e-7);
}

NOVA_TEST(elasticity_3d_rigid_translation_and_rotation_have_zero_strain) {
    const auto box = make_box_mesh_3d(
        0.0,1.0,0.0,1.0,0.0,1.0,2U,2U,2U);
    auto problem = make_model(box, 2000.0, 0.25);
    const auto rigid = [](const Point3D& p) -> ElasticVector3D {
        return {0.4 - 0.02*p.y, -0.3 + 0.02*p.x, 0.2};
    };
    fix_exterior(problem, box.boundary_faces, rigid);
    const auto solved = solve_linear_elasticity_3d(problem);
    NOVA_REQUIRE(solved.linear_solver.converged);
    for (std::size_t i = 0; i < solved.displacement.size(); ++i) {
        const auto expected = rigid(problem.mesh.nodes[i]);
        for (std::size_t c = 0; c < 3U; ++c) {
            NOVA_REQUIRE_NEAR(solved.displacement[i][c], expected[c], 1.0e-9);
        }
    }
    for (const auto& strain : solved.cell_strain) {
        for (const double component : strain) {
            NOVA_REQUIRE(std::abs(component) < 1.0e-8);
        }
    }
    NOVA_REQUIRE(std::abs(solved.strain_energy_joules) < 1.0e-12);
}

NOVA_TEST(elasticity_3d_uniaxial_traction_recovers_analytic_stress) {
    const auto box = make_box_mesh_3d(
        0.0,1.0,0.0,1.0,0.0,1.0,2U,2U,2U);
    auto problem = make_model(box, 100.0, 0.0);
    for (std::size_t node = 0; node < box.mesh.nodes.size(); ++node) {
        if (box.mesh.nodes[node].x == 0.0) {
            problem.prescribed_displacements.push_back({node,{0.0,0.0,0.0}});
        }
    }
    for (const auto& face : box.boundary_faces) {
        if (box.mesh.nodes[face[0]].x == 1.0 &&
            box.mesh.nodes[face[1]].x == 1.0 &&
            box.mesh.nodes[face[2]].x == 1.0) {
            problem.face_tractions.push_back({face,{10.0,0.0,0.0}});
        }
    }
    NOVA_REQUIRE(!problem.face_tractions.empty());
    const auto solved = solve_linear_elasticity_3d(problem);
    NOVA_REQUIRE(solved.linear_solver.converged);
    for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
        const auto& p = problem.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(solved.displacement[i][0], 0.1*p.x, 1.0e-8);
        NOVA_REQUIRE_NEAR(solved.displacement[i][1], 0.0, 1.0e-8);
        NOVA_REQUIRE_NEAR(solved.displacement[i][2], 0.0, 1.0e-8);
    }
    for (const auto& stress : solved.cell_stress) {
        NOVA_REQUIRE_NEAR(stress[0], 10.0, 1.0e-6);
        for (std::size_t c = 1; c < 6U; ++c) {
            NOVA_REQUIRE_NEAR(stress[c], 0.0, 1.0e-6);
        }
    }
    NOVA_REQUIRE_NEAR(solved.strain_energy_joules, 0.5, 1.0e-7);
}

NOVA_TEST(elasticity_3d_orientation_independent_constitutive_recovery) {
    LinearElasticity3DProblem problem;
    problem.mesh.nodes = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
    problem.mesh.elements = {{0,1,2,3}};
    problem.cell_materials = {{2000.0,0.25}};
    for (std::size_t node=0; node<4U; ++node) {
        problem.prescribed_displacements.push_back({
            node, affine_displacement(problem.mesh.nodes[node])});
    }
    const auto forward = solve_linear_elasticity_3d(problem);
    std::swap(problem.mesh.elements[0][1], problem.mesh.elements[0][2]);
    const auto reversed = solve_linear_elasticity_3d(problem);
    NOVA_REQUIRE(forward.linear_solver.converged);
    NOVA_REQUIRE(reversed.linear_solver.converged);
    NOVA_REQUIRE_NEAR(forward.strain_energy_joules,
                      reversed.strain_energy_joules, 1.0e-12);
    for (std::size_t i=0; i<6U; ++i) {
        NOVA_REQUIRE_NEAR(forward.cell_strain[0][i],
                          reversed.cell_strain[0][i], 1.0e-12);
        NOVA_REQUIRE_NEAR(forward.cell_stress[0][i],
                          reversed.cell_stress[0][i], 1.0e-9);
    }
}

NOVA_TEST(elasticity_3d_rejects_missing_materials_bad_anchors_and_faces) {
    const auto box = make_box_mesh_3d(
        0.0,1.0,0.0,1.0,0.0,1.0,2U,2U,2U);
    auto problem = make_model(box, 100.0, 0.25);
    fix_exterior(problem, box.boundary_faces, affine_displacement);

    auto bad = problem;
    bad.cell_materials.clear();
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(bad); }));
    bad = problem;
    bad.cell_materials[0].poisson_ratio = 0.5;
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(bad); }));
    bad = problem;
    bad.cell_materials[0].young_modulus = -1.0;
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(bad); }));
    bad = problem;
    bad.prescribed_displacements.resize(2);
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(bad); }));
    bad = problem;
    bad.prescribed_displacements = {
        {0U,{0,0,0}}, {1U,{0,0,0}}, {2U,{0,0,0}}};
    // Independent of the box numbering, three consecutive nodes might
    // not be collinear; enforce a known line x=0, y=0 across the z axis.
    bad.prescribed_displacements.clear();
    for (std::size_t node=0; node<bad.mesh.nodes.size(); ++node) {
        const auto& p=bad.mesh.nodes[node];
        if (p.x==0.0 && p.y==0.0) {
            bad.prescribed_displacements.push_back({node,{0,0,0}});
        }
    }
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(bad); }));
    bad = problem;
    bad.face_tractions.push_back({box.boundary_faces.front(), {1,0,0}});
    bad.face_tractions.push_back({box.boundary_faces.front(), {1,0,0}});
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(bad); }));
    bad = problem;
    bad.cell_body_force_density.assign(bad.mesh.elements.size(), {0,0,0});
    bad.cell_body_force_density[0][0] =
        std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(bad); }));
}


NOVA_TEST(elasticity_3d_uniform_thermal_expansion_is_stress_free) {
    const auto box = make_box_mesh_3d(
        0.0,1.0,0.0,1.0,0.0,1.0,2U,2U,2U);
    auto problem = make_model(box, 1200.0, 0.25);
    problem.cell_thermal_expansion_per_kelvin.assign(
        problem.mesh.elements.size(), 0.002);
    problem.cell_temperature_change_kelvin.assign(
        problem.mesh.elements.size(), 5.0); // alpha * deltaT = 0.01
    const auto free_dilation = [](const Point3D& p) -> ElasticVector3D {
        return {0.01*p.x, 0.01*p.y, 0.01*p.z};
    };
    fix_exterior(problem, box.boundary_faces, free_dilation);
    const auto solved = solve_linear_elasticity_3d(problem);
    NOVA_REQUIRE(solved.linear_solver.converged);
    for (std::size_t n=0; n<problem.mesh.nodes.size(); ++n) {
        const auto expected = free_dilation(problem.mesh.nodes[n]);
        for (std::size_t c=0; c<3U; ++c) {
            NOVA_REQUIRE_NEAR(solved.displacement[n][c], expected[c], 1.0e-8);
        }
    }
    for (std::size_t e=0; e<solved.cell_stress.size(); ++e) {
        const auto& total_strain = solved.cell_strain[e];
        for (std::size_t c=0; c<3U; ++c) {
            NOVA_REQUIRE_NEAR(total_strain[c], 0.01, 1.0e-8);
        }
        for (std::size_t c=3; c<6U; ++c) {
            NOVA_REQUIRE_NEAR(total_strain[c], 0.0, 1.0e-8);
        }
        for (const double stress : solved.cell_stress[e]) {
            NOVA_REQUIRE_NEAR(stress, 0.0, 1.0e-5);
        }
    }
    NOVA_REQUIRE_NEAR(solved.strain_energy_joules, 0.0, 1.0e-8);
}

NOVA_TEST(elasticity_3d_fully_clamped_thermal_stress_and_stored_energy) {
    const auto box = make_box_mesh_3d(
        0.0,1.0,0.0,1.0,0.0,1.0,2U,2U,2U);
    auto problem = make_model(box, 1200.0, 0.25);
    problem.cell_thermal_expansion_per_kelvin.assign(
        problem.mesh.elements.size(), 0.002);
    problem.cell_temperature_change_kelvin.assign(
        problem.mesh.elements.size(), 5.0);
    // E=1200, nu=0.25 => lambda=mu=480, 3*lambda+2*mu=2400.
    // Clamped alpha*deltaT=0.01 gives sigma_ii=-24 Pa and U=0.36 J.
    for (std::size_t node=0; node<problem.mesh.nodes.size(); ++node) {
        problem.prescribed_displacements.push_back({node,{0.0,0.0,0.0}});
    }
    const auto solved = solve_linear_elasticity_3d(problem);
    NOVA_REQUIRE(solved.linear_solver.converged);
    for (const auto& u : solved.displacement) {
        for (const double component : u) {
            NOVA_REQUIRE_NEAR(component, 0.0, 1.0e-12);
        }
    }
    for (std::size_t e=0; e<solved.cell_stress.size(); ++e) {
        for (std::size_t c=0; c<6U; ++c) {
            NOVA_REQUIRE_NEAR(solved.cell_strain[e][c], 0.0, 1.0e-12);
            NOVA_REQUIRE_NEAR(solved.cell_stress[e][c],
                              c<3U ? -24.0 : 0.0, 1.0e-8);
        }
    }
    NOVA_REQUIRE_NEAR(solved.strain_energy_joules, 0.36, 1.0e-9);
}

NOVA_TEST(elasticity_3d_rejects_malformed_or_nonfinite_thermal_inputs) {
    const auto box = make_box_mesh_3d(
        0.0,1.0,0.0,1.0,0.0,1.0,1U,1U,1U);
    auto problem = make_model(box, 1200.0, 0.25);
    const auto zero = [](const Point3D&) -> ElasticVector3D {
        return {0.0,0.0,0.0};
    };
    fix_exterior(problem, box.boundary_faces, zero);
    const auto cell_count = problem.mesh.elements.size();
    problem.cell_thermal_expansion_per_kelvin.assign(cell_count, 0.002);
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(problem); }));
    problem.cell_temperature_change_kelvin.assign(cell_count - 1U, 5.0);
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(problem); }));
    problem.cell_temperature_change_kelvin.assign(cell_count, 5.0);
    problem.cell_thermal_expansion_per_kelvin[0] =
        std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(problem); }));
    problem.cell_thermal_expansion_per_kelvin[0] =
        std::numeric_limits<double>::max();
    NOVA_REQUIRE(rejects([&] { return solve_linear_elasticity_3d(problem); }));
    // Zero thermal strain must exactly preserve the source-free path.
    problem.cell_thermal_expansion_per_kelvin.assign(cell_count, 0.0);
    const auto solved = solve_linear_elasticity_3d(problem);
    NOVA_REQUIRE(solved.linear_solver.converged);
    NOVA_REQUIRE_NEAR(solved.strain_energy_joules, 0.0, 1.0e-12);
}
