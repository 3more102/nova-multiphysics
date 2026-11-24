#include "test_framework.hpp"
#include "nova/physics/maxwell_power_audit_3d.hpp"
#include <complex>
#include <stdexcept>

using namespace nova::physics;

NOVA_TEST(maxwell_power_audit_lossy_manufactured_field) {
    using C=std::complex<double>;
    const auto mesh=nova::fem::make_box_mesh_3d(0,1,0,1,0,1,1,1,1).mesh;
    constexpr double k=2.3, c=299792458.0, eps0=8.8541878128e-12;
    const ComplexVector3D exact{C(1,0.2),C(-0.3,0.4),C(0.5,-0.1)};
    std::vector<MaxwellCellMaterial3D> materials;
    std::vector<ComplexVector3D> sources;
    double mean_loss=0.0;
    for (std::size_t i=0;i<mesh.elements.size();++i) {
        const C eps(i%2 ? -2.0 : 3.0,i%2 ? 0.7 : 0.2);
        materials.push_back({eps,1.0});
        ComplexVector3D f{};
        for (std::size_t d=0;d<3;++d) f[d]=-k*k*eps*exact[d];
        sources.push_back(f);
        mean_loss+=eps.imag()/mesh.elements.size();
    }
    const auto system=assemble_maxwell_3d(mesh,materials,k,false);
    const auto b=assemble_maxwell_cell_load_3d(mesh,sources,system);
    nova::linalg::GmresOptions opts;
    opts.restart=100;
    opts.solver.max_iterations=1000;
    opts.solver.relative_tolerance=1e-12;
    const auto sol=solve_maxwell_3d(system,b,opts);
    NOVA_REQUIRE(sol.diagnostics.converged);
    const auto a=audit_maxwell_power_3d(mesh,materials,k,system,b,sol.edge_field);
    const double expected=0.5*k*c*eps0*mean_loss*
        (std::norm(exact[0])+std::norm(exact[1])+std::norm(exact[2]));
    NOVA_REQUIRE_NEAR(a.source_power,expected,1e-10);
    NOVA_REQUIRE_NEAR(a.operator_absorbed_power,expected,1e-10);
    NOVA_REQUIRE_NEAR(a.cell_absorbed_power,expected,1e-10);
    NOVA_REQUIRE(a.relative_power_imbalance<1e-8);
    NOVA_REQUIRE(a.relative_operator_cell_mismatch<1e-8);
    NOVA_REQUIRE_NEAR(a.residual_work,0.0,1e-10);
    auto altered=b;
    for (auto& z : altered) z*=2.0;
    const auto bad=audit_maxwell_power_3d(mesh,materials,k,system,altered,sol.edge_field);
    NOVA_REQUIRE_NEAR(bad.source_power,2*expected,1e-10);
    NOVA_REQUIRE_NEAR(bad.relative_power_imbalance,0.5,1e-8);
    NOVA_REQUIRE_NEAR(bad.residual_work,expected,1e-10);
    NOVA_REQUIRE(bad.algebraic_residual_norm>1e-3);
    NOVA_REQUIRE(bad.relative_operator_cell_mismatch<1e-8);
}

NOVA_TEST(maxwell_power_audit_zero_and_pec_validation) {
    const auto mesh=nova::fem::make_box_mesh_3d(0,1,0,1,0,1,1,1,1).mesh;
    const std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size(),{{2,0},1});
    const auto s=assemble_maxwell_3d(mesh,materials,2.3,true);
    std::vector<std::complex<double>> b(s.edges.size()), x(s.edges.size());
    const auto a=audit_maxwell_power_3d(mesh,materials,2.3,s,b,x);
    NOVA_REQUIRE_NEAR(a.relative_power_imbalance,0.0,0.0);
    NOVA_REQUIRE_NEAR(a.relative_operator_cell_mismatch,0.0,0.0);
    NOVA_REQUIRE_NEAR(a.algebraic_residual_norm,0.0,0.0);
    for (std::size_t i=0;i<x.size();++i)
        if (i!=s.free_edges.front()) { x[i]=1.0; break; }
    bool rejected=false;
    try { (void)audit_maxwell_power_3d(mesh,materials,2.3,s,b,x); }
    catch(const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected);
}
