#include "test_framework.hpp"

#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/error_norms_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"
#include "nova/io/vtk_legacy_3d.hpp"
#include "nova/physics/thermoelasticity_3d.hpp"

#include <algorithm>
#include <array>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>

namespace {
using namespace nova::fem;
using namespace nova::physics;

template<class Action>
bool rejects(Action action) {
    try { action(); } catch (const std::exception&) { return true; }
    return false;
}

TetrahedralMesh3D unit_tet() {
    return {{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}, {{0,1,2,3}}};
}
const std::array<DisplacementConstraint3D, 3> base{{{0,{}},{1,{}},{2,{}}}};
const std::array<IsotropicElasticMaterial3D, 1> elastic{{{100,0}}};
const std::array<ThermalExpansion3D, 1> expansion{{{0.01,300}}};

} // namespace

NOVA_TEST(thermoelastic_3d_p1_trace_load_stress_and_exact_quadratic_energy) {
    auto mesh = unit_tet();
    for (unsigned orientation = 0; orientation < 2; ++orientation) {
        PreparedElasticity3D op(mesh, elastic, base);
        const std::array<double, 4> temperature{300,300,300,310};
        const auto r = op.solve_temperature(temperature, expansion);
        NOVA_REQUIRE(r.linear_solver.converged);
        NOVA_REQUIRE_NEAR(r.displacement[3][2], 0.025, 1e-13);
        const ElasticTensor6 stress{-2.5,-2.5,0,0,0,0};
        for (std::size_t i = 0; i < 6; ++i) {
            NOVA_REQUIRE_NEAR(r.cell_stress[0][i], stress[i], 1e-12);
        }
        // e_zz=1/40, s(z)=z/10; integrate 50*(2*s^2+(1/40-s)^2).
        // int(1,z,z^2)=(1/6,1/24,1/60), independent of the transfer rule.
        constexpr double exact = 50.0*(0.03/60.0 - 0.005/24.0 + 0.000625/6.0);
        NOVA_REQUIRE_NEAR(r.strain_energy_joules, exact, 1e-14);
        const std::array<double, 1> mean{0.025};
        const auto p0 = op.solve_cell_eigenstrain(mean);
        NOVA_REQUIRE_NEAR(p0.displacement[3][2], 0.025, 1e-13);
        NOVA_REQUIRE_NEAR(p0.strain_energy_joules, 0.0625/6.0, 1e-14);
        NOVA_REQUIRE(r.strain_energy_joules > p0.strain_energy_joules);
        std::swap(mesh.elements[0][0], mesh.elements[0][1]);
    }
}

NOVA_TEST(thermoelastic_3d_owned_snapshot_ic0_reuse_and_failure_recovery) {
    auto box = make_box_mesh_3d(0,1,0,1,0,1,2,2,2);
    auto mats = std::vector<IsotropicElasticMaterial3D>(box.mesh.elements.size(), {100,0.25});
    const std::vector<ThermalExpansion3D> laws(mats.size(), {0.001,300});
    std::vector<DisplacementConstraint3D> fixed;
    for (std::size_t i = 0; i < box.mesh.nodes.size(); ++i) {
        if (box.mesh.nodes[i].x == 0) fixed.push_back({i,{0,0,0}});
    }
    std::vector<double> temp(box.mesh.nodes.size(), 305);
    PreparedElasticity3D op(box.mesh, mats, fixed);
    box = {};
    mats.clear();
    nova::linalg::SolverOptions opts;
    opts.relative_tolerance = 1e-12;
    opts.preconditioner = nova::linalg::PreconditionerType::IncompleteCholesky0;
    const auto heated = op.solve_temperature(temp, laws, opts);
    NOVA_REQUIRE(heated.linear_solver.converged && op.factor_builds() == 1);
    NOVA_REQUIRE(op.stiffness_nonzeros() > 3U*temp.size());
    std::fill(temp.begin(), temp.end(), 295);
    auto limited = opts;
    limited.max_iterations = 1;
    NOVA_REQUIRE(!op.solve_temperature(temp, laws, limited).linear_solver.converged);
    const auto cooled = op.solve_temperature(temp, laws, opts);
    NOVA_REQUIRE(cooled.linear_solver.converged && op.factor_builds() == 1);
    for (std::size_t i = 0; i < temp.size(); ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            NOVA_REQUIRE_NEAR(cooled.displacement[i][j], -heated.displacement[i][j], 2e-12);
        }
    }
    NOVA_REQUIRE_NEAR(cooled.strain_energy_joules, heated.strain_energy_joules, 2e-12);
}

NOVA_TEST(thermoelastic_3d_rejects_invalid_fields_constraints_and_geometry) {
    auto mesh = unit_tet();
    PreparedElasticity3D op(mesh, elastic, base);
    std::vector<double> field(4, 300);
    NOVA_REQUIRE(rejects([&] { (void)op.solve_temperature({}, expansion); }));
    NOVA_REQUIRE(rejects([&] { (void)op.solve_temperature(field, {}); }));
    auto bad = expansion;
    bad.front().reference_temperature = std::numeric_limits<double>::infinity();
    NOVA_REQUIRE(rejects([&] { (void)op.solve_temperature(field, bad); }));
    field[2] = std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects([&] { (void)op.solve_temperature(field, expansion); }));
    NOVA_REQUIRE(rejects([&] { PreparedElasticity3D invalid(mesh, elastic, {}); }));
    mesh.nodes.back() = {0.2,0.2,0};
    NOVA_REQUIRE(rejects([&] { PreparedElasticity3D invalid(mesh, elastic, base); }));
}

NOVA_TEST(thermoelastic_3d_requires_absolute_steady_temperature_on_every_component) {
    ScalarDiffusion3DProblem p;
    p.mesh = unit_tet();
    p.mesh.nodes.insert(p.mesh.nodes.end(), {{2,0,0},{3,0,0},{2,1,0},{2,0,1}});
    p.mesh.elements.push_back({4,5,6,7});
    p.cell_diffusivity = {2, 2};
    p.dirichlet_constraints.push_back({0,300});
    NOVA_REQUIRE(solve_scalar_diffusion_3d(p).linear_solver.converged);
    p.require_anchored_components = true;
    NOVA_REQUIRE(rejects([&] { (void)solve_scalar_diffusion_3d(p); }));
    p.robin_faces.push_back({{4,5,6}, [](double,double,double){return 0.;},
                                       [](double,double,double){return 270.;}});
    NOVA_REQUIRE(rejects([&] { (void)solve_scalar_diffusion_3d(p); }));
    p.robin_faces[0].transfer_coefficient = [](double,double,double){return 3.;};
    const auto result = solve_scalar_diffusion_3d(p);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t i = 0; i < result.values.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.values[i], i < 4 ? 300 : 270, 1e-8);
    }
}

NOVA_TEST(thermoelastic_3d_manufactured_quadratic_displacement_convergence) {
    // T=300+x, alpha=1. u=(x^2-y^2-z^2,2xy,2xz)/2 has eps=x I
    // and zero exact stress. It admits nontrivial unconstrained interior DOFs.
    const auto x = [](double a,double,double){return a;};
    const auto y = [](double,double b,double){return b;};
    const auto z = [](double,double,double c){return c;};
    const auto zero = [](double,double,double){return 0.;};
    const std::array<ExactScalarField3D,3> exact{{
        {[](double a,double b,double c){return 0.5*(a*a-b*b-c*c);}, x,
         [](double,double b,double){return -b;}, [](double,double,double c){return -c;}},
        {[](double a,double b,double){return a*b;}, y, x, zero},
        {[](double a,double,double c){return a*c;}, z, zero, x}
    }};
    double prior_l2 = 0, prior_h1 = 0, prior_energy = 0;
    for (const std::size_t n : {2U,4U,8U}) {
        auto box = make_box_mesh_3d(0,1,0,1,0,1,n,n,n);
        ScalarDiffusion3DProblem heat;
        heat.mesh = std::move(box.mesh);
        heat.diffusivity_tensor = [](double,double,double){return SymmetricTensor3D{2,3,4,0.3,-0.4,0.2};};
        std::set<std::size_t> boundary;
        for (const auto& face : box.boundary_faces) {
            heat.dirichlet_faces.push_back({face, [](double a,double,double){return 300.+a;}});
            boundary.insert(face.begin(),face.end());
        }
        const auto temperature = solve_scalar_diffusion_3d(heat, {1000,1e-13,1e-14});
        NOVA_REQUIRE(temperature.linear_solver.converged);
        const auto field = evaluate_p1_tetrahedral_flux_3d(heat, temperature.values);
        for (const auto& cell : field) {
            NOVA_REQUIRE_NEAR(cell.gradient.x, 1, 3e-9);
            NOVA_REQUIRE_NEAR(cell.flux.x, -2, 8e-9);
            NOVA_REQUIRE_NEAR(cell.flux.y, -0.3, 8e-9);
            NOVA_REQUIRE_NEAR(cell.flux.z, 0.4, 8e-9);
        }
        std::vector<DisplacementConstraint3D> fixed;
        for (const auto node : boundary) {
            const auto [a,b,c] = heat.mesh.nodes[node];
            fixed.push_back({node, {exact[0].value(a,b,c),exact[1].value(a,b,c),exact[2].value(a,b,c)}});
        }
        const std::vector<IsotropicElasticMaterial3D> material(heat.mesh.elements.size(), {10,0.25});
        const std::vector<ThermalExpansion3D> laws(material.size(), {1,300});
        PreparedElasticity3D op(heat.mesh, material, fixed);
        const auto mech = op.solve_temperature(temperature.values, laws,
            {1000,1e-12,1e-14,nova::linalg::PreconditionerType::Jacobi});
        NOVA_REQUIRE(mech.linear_solver.converged);
        double l2 = 0, h1 = 0;
        std::vector<double> component(heat.mesh.nodes.size());
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < component.size(); ++j) component[j] = mech.displacement[j][i];
            const auto norms = compute_error_norms_3d(heat.mesh, component, exact[i]);
            l2 = std::hypot(l2, norms.l2);
            h1 = std::hypot(h1, norms.h1_seminorm);
        }
        const double energy = mech.strain_energy_joules;
        std::cout << std::setprecision(10) << "coupled refinement n=" << n << " L2_u=" << l2
                  << " H1semi_u=" << h1 << " energy_J=" << energy << '\n';
        if (prior_l2 > 0) {
            // Factor-two refinement: expected L2=O(h^2), H1=O(h), energy=O(h^2).
            NOVA_REQUIRE(std::log2(prior_l2/l2) > 1.8);
            NOVA_REQUIRE(std::log2(prior_h1/h1) > 0.9);
            NOVA_REQUIRE(std::log2(prior_energy/energy) > 1.8);
        }
        prior_l2 = l2;
        prior_h1 = h1;
        prior_energy = energy;
    }
}

NOVA_TEST(thermoelastic_vtk_rejects_bad_fields_before_writing) {
    auto mesh = unit_tet();
    PreparedElasticity3D op(mesh, elastic, base);
    const std::vector<double> t{300,300,300,310};
    const std::array<int,1> tags{42};
    std::array<P1TetrahedralFlux3D,1> fields{};
    auto state = op.solve_temperature(t, expansion);
    std::ostringstream out;
    const auto write = [&] { nova::io::write_legacy_vtk_thermoelastic_3d(out,mesh,t,tags,fields,state); };
    state.linear_solver.converged = false;
    NOVA_REQUIRE(rejects(write) && out.str().empty());
    state.linear_solver.converged = true;
    fields[0].gradient.x = std::numeric_limits<double>::infinity();
    NOVA_REQUIRE(rejects(write) && out.str().empty());
    fields[0].gradient.x = 0;
    state.cell_stress.clear();
    NOVA_REQUIRE(rejects(write) && out.str().empty());
}

NOVA_TEST(thermoelastic_3d_partial_constraints_preserve_free_expansion_and_reject_a_hinge) {
    auto mesh = unit_tet();
    std::vector<DisplacementComponent3D> supports{{0,0,0.2},{0,1,0},{0,2,0},
                                                   {1,1,0},{1,2,0},{2,2,0}};
    PreparedElasticity3D free(mesh,elastic,{}, {}, {},supports);
    const std::array<double,4> t{305,305,305,305};
    const auto result = free.solve_temperature(t,expansion);
    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t n = 0; n < t.size(); ++n) {
        const auto [x,y,z] = mesh.nodes[n];
        const std::array<double,3> exact{0.2+.05*x,.05*y,.05*z};
        for (std::size_t i = 0; i < 3; ++i) NOVA_REQUIRE_NEAR(result.displacement[n][i],exact[i],1e-12);
    }
    NOVA_REQUIRE_NEAR(result.strain_energy_joules,0,1e-20);
    supports.pop_back(); // Free rotation about x.
    NOVA_REQUIRE(rejects([&] { PreparedElasticity3D loose(mesh,elastic,{}, {}, {},supports); }));
    // A second tet shares only the x edge; it can still rotate about it.
    mesh.nodes.insert(mesh.nodes.end(), {{0,-1,0},{0,0,-1}});
    mesh.elements.push_back({0,1,4,5});
    const auto topology = analyze_tetrahedral_topology_3d(mesh);
    NOVA_REQUIRE(topology.cell_component == std::vector<std::size_t>({0,1}));
    const std::array<IsotropicElasticMaterial3D,2> materials{elastic[0],elastic[0]};
    NOVA_REQUIRE(rejects([&] { PreparedElasticity3D hinge(mesh,materials,base); }));
}

NOVA_TEST(thermoelastic_vtk_serializes_tensor_shear_in_cartesian_order) {
    const auto mesh = unit_tet();
    const std::vector<double> t(4,300);
    std::vector<DisplacementConstraint3D> prescribed;
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
        const auto& p = mesh.nodes[i];
        prescribed.push_back({i,{.02*p.y,.04*p.z,.06*p.x}});
    }
    PreparedElasticity3D op(mesh,elastic,prescribed);
    const auto r = op.solve_temperature(t,expansion);
    std::ostringstream out;
    const std::array<int,1> tags{42};
    const std::array<P1TetrahedralFlux3D,1> fields{};
    nova::io::write_legacy_vtk_thermoelastic_3d(out,mesh,t,tags,fields,r);
    std::istringstream data(out.str());
    std::string line;
    while (std::getline(data,line) && line != "TENSORS mean_stress double") {}
    NOVA_REQUIRE(line == "TENSORS mean_stress double");
    // E=100, nu=0: sigma_xy=1, sigma_yz=2, sigma_xz=3 Pa.
    for (const double expected : {0.,1.,3.,1.,0.,2.,3.,2.,0.}) {
        double value{};
        NOVA_REQUIRE(static_cast<bool>(data >> value));
        NOVA_REQUIRE_NEAR(value,expected,1e-12);
    }
}

NOVA_TEST(thermoelastic_3d_small_load_change_respects_physical_force_tolerance) {
    const auto mesh = make_box_mesh_3d(0,1,0,1,0,1,3,3,3).mesh;
    const std::vector<IsotropicElasticMaterial3D> material(mesh.elements.size(), {1e11,.29});
    const std::vector<ThermalExpansion3D> laws(material.size(), {1.2e-5,300});
    std::vector<DisplacementConstraint3D> fixed;
    for (std::size_t n = 0; n < mesh.nodes.size(); ++n) {
        if (mesh.nodes[n].x == 0) fixed.push_back({n,{}});
    }
    PreparedElasticity3D op(mesh,material,fixed);
    std::vector<double> temp(mesh.nodes.size(),305);
    nova::linalg::SolverOptions opts{1000,1e-11,0,nova::linalg::PreconditionerType::Jacobi};
    const auto loaded = op.solve_temperature(temp,laws,opts);
    NOVA_REQUIRE(loaded.linear_solver.converged);
    const auto repeated = op.solve_temperature(temp,laws,opts);
    NOVA_REQUIRE(repeated.linear_solver.converged && repeated.linear_solver.iterations == 0);
    std::fill(temp.begin(),temp.end(),300);
    const auto unloaded = op.solve_temperature(temp,laws,opts);
    NOVA_REQUIRE(unloaded.linear_solver.converged);
    NOVA_REQUIRE_NEAR(unloaded.strain_energy_joules,0,1e-24);
}
