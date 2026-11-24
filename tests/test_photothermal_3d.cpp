#include "test_framework.hpp"
#include "nova/physics/photothermal_3d.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
using namespace nova::fem;
using namespace nova::physics;
Photothermal3DProblem plate(double loss, double amplitude) {
    const auto box = make_box_mesh_3d(0,1,0,1,0,1,2,2,2);
    Photothermal3DProblem p;
    p.thermal.mesh = box.mesh;
    p.thermal.cell_diffusivity.assign(box.mesh.elements.size(), 4.0);
    p.vacuum_wavenumber = 2.3;
    p.pec_on_exterior = false;
    const std::complex<double> eps(2.0, loss);
    p.optical_materials.assign(box.mesh.elements.size(), {eps,1.0});
    ComplexVector3D source{};
    source[0] = -p.vacuum_wavenumber*p.vacuum_wavenumber*eps*amplitude;
    p.optical_cell_sources.assign(box.mesh.elements.size(), source);
    for (const auto& face : box.boundary_faces) {
        const auto on_x = [&](double x) {
            return std::all_of(face.begin(), face.end(),
                [&](std::size_t node) { return box.mesh.nodes[node].x == x; });
        };
        if (on_x(0.0) || on_x(1.0))
            p.thermal.dirichlet_faces.push_back(
                {face, [](double,double,double) { return 300.0; }});
    }
    return p;
}

nova::linalg::GmresOptions accurate_optical_options() {
    nova::linalg::GmresOptions o;
    o.restart=100;
    o.solver.max_iterations=1000;
    o.solver.relative_tolerance=1e-12;
    return o;
}
}

NOVA_TEST(photothermal_constant_field_power_and_thermal_balance) {
    auto p=plate(0.6,100.0);
    p.thermal.cell_source.assign(p.thermal.mesh.elements.size(), 7.0);
    p.thermal.source=[](double,double,double){return 3.0;};
    constexpr double c=299792458.0, eps0=8.8541878128e-12;
    const double expected=0.5*p.vacuum_wavenumber*c*eps0*0.6*10000.0;
    const auto r=solve_one_way_photothermal_3d(p,accurate_optical_options());
    NOVA_REQUIRE(r.optical_solution.diagnostics.converged);
    NOVA_REQUIRE(r.thermal_solution.linear_solver.converged);
    NOVA_REQUIRE_NEAR(r.total_absorbed_power_watts,expected,1e-7);
    NOVA_REQUIRE_NEAR(r.thermal_balance.generated_heat,expected+10.0,1e-7);
    NOVA_REQUIRE_NEAR(r.thermal_balance.signed_balance_defect,0.0,1e-6);
    NOVA_REQUIRE_NEAR(r.thermal_balance.constraint_outward_heat,expected+10.0,1e-6);
    for (double q : r.absorbed_power_density) NOVA_REQUIRE_NEAR(q,expected,1e-7);
    for (double q : p.thermal.cell_source) NOVA_REQUIRE_NEAR(q,7.0,0.0);
    NOVA_REQUIRE(*std::max_element(r.thermal_solution.values.begin(),
        r.thermal_solution.values.end())>300.0);
}

NOVA_TEST(photothermal_lossless_zero_absorption) {
    const auto p=plate(0.0,0.0);
    const auto r=solve_one_way_photothermal_3d(p,accurate_optical_options());
    NOVA_REQUIRE_NEAR(r.total_absorbed_power_watts,0.0,0.0);
    for (double t : r.thermal_solution.values) NOVA_REQUIRE_NEAR(t,300.0,1e-7);
}

NOVA_TEST(photothermal_rejects_mismatched_cell_data) {
    auto p=plate(0.6,100.0);
    p.optical_cell_sources.pop_back();
    bool rejected=false;
    try { (void)solve_one_way_photothermal_3d(p); }
    catch (const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected);
}
