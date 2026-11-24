#include "nova/physics/optical_helmholtz_2d.hpp"
#include "test_framework.hpp"
#include <algorithm>
#include <complex>
#include <numeric>
#include <stdexcept>

using namespace nova::optics;
using C=std::complex<double>;
namespace {
double wave_error(std::size_t nx, ScalarPolarization2D pol, C eps, C mu, bool reverse) {
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,0.25,nx,nx/4);
    if (reverse) for (auto& e:rect.mesh.elements) std::swap(e[1],e[2]);
    const double k=2.0;
    const C a=1.0/(pol==ScalarPolarization2D::TE_Ez?mu:eps);
    const C q=k*std::sqrt(eps*mu), I(0,1);
    std::vector<HelmholtzMaterial2D> mat(rect.mesh.elements.size(),{eps,mu});
    std::vector<HelmholtzRobin2D> bc;
    for (auto e:rect.boundary_edges) {
        const auto x=rect.mesh.nodes[e[0]].x, y=rect.mesh.nodes[e[1]].x;
        if (x==0 && y==0) bc.push_back({e,-I*a*q,-2.0*I*a*q});
        if (x==1 && y==1) bc.push_back({e,-I*a*q,0});
    }
    const auto result=solve_helmholtz_2d(rect.mesh,k,pol,mat,{},bc);
    double error=0;
    for (std::size_t i=0;i<result.field.size();++i)
        error=std::max(error,std::abs(result.field[i]-std::exp(I*q*rect.mesh.nodes[i].x)));
    NOVA_REQUIRE(result.relative_backward_error<1e-12);
    const double loss=std::accumulate(result.cell_loss.begin(),result.cell_loss.end(),0.0);
    NOVA_REQUIRE_NEAR(result.material_loss,loss,1e-13);
    NOVA_REQUIRE(result.excitation_input>0);
    NOVA_REQUIRE(result.robin_outflow>0);
    NOVA_REQUIRE_NEAR(result.dirichlet_input,0,1e-14);
    NOVA_REQUIRE_NEAR(result.power_balance_residual,0,1e-10);
    const C b=pol==ScalarPolarization2D::TE_Ez?eps:mu;
    const double decay=2*q.imag();
    const double integral=decay==0?0.25:0.25*(-std::expm1(-decay))/decay;
    const double expected=(-a.imag()*std::norm(q)+k*k*b.imag())*integral;
    if (expected>0) NOVA_REQUIRE(std::abs(loss-expected)/expected<0.04);
    else NOVA_REQUIRE_NEAR(loss,0,1e-14);
    return error;
}
}
NOVA_TEST(helmholtz_2d_plane_wave_te_tm_convergence) {
    for (auto p:{ScalarPolarization2D::TE_Ez,ScalarPolarization2D::TM_Hz}) {
        const double coarse=wave_error(12,p,2.25,1.0,false);
        const double fine=wave_error(24,p,2.25,1.0,false);
        NOVA_REQUIRE(fine<0.01);
        NOVA_REQUIRE(fine<0.4*coarse);
        NOVA_REQUIRE_NEAR(wave_error(24,p,2.25,1.0,true),fine,1e-12);
    }
}
NOVA_TEST(helmholtz_2d_lossy_wave_electric_and_magnetic) {
    for (auto p:{ScalarPolarization2D::TE_Ez,ScalarPolarization2D::TM_Hz})
        NOVA_REQUIRE(wave_error(32,p,C(2.25,0.2),C(1.0,0.1),false)<0.01);
}
NOVA_TEST(helmholtz_2d_enz_te_affine_tm_rejected) {
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,1,4,3);
    std::vector<HelmholtzMaterial2D> mat(rect.mesh.elements.size(),{0.0,1.0});
    std::vector<HelmholtzDirichlet2D> bc;
    for (auto i:rect.boundary_nodes) {
        auto x=rect.mesh.nodes[i].x;
        if (x==0 || x==1) bc.push_back({i,C(1+x,2*x)});
    }
    const auto result=solve_helmholtz_2d(rect.mesh,2,ScalarPolarization2D::TE_Ez,mat,bc);
    for (std::size_t i=0;i<result.field.size();++i) {
        const auto x=rect.mesh.nodes[i].x;
        NOVA_REQUIRE_NEAR(std::abs(result.field[i]-C(1+x,2*x)),0,1e-12);
    }
    bool rejected=false;
    try { (void)solve_helmholtz_2d(rect.mesh,2,ScalarPolarization2D::TM_Hz,mat,bc); }
    catch (const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected);
}
NOVA_TEST(helmholtz_2d_rejects_interior_robin_and_active_material) {
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,1,1,1);
    std::vector<HelmholtzMaterial2D> mat(rect.mesh.elements.size());
    std::vector<HelmholtzRobin2D> bc{{{0,3},0,1}};
    bool rejected=false;
    try { (void)solve_helmholtz_2d(rect.mesh,2,ScalarPolarization2D::TE_Ez,mat,{},bc); }
    catch (const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected);
    mat[0].epsilon_r=C(1,-0.1);
    rejected=false;
    try { (void)solve_helmholtz_2d(rect.mesh,2,ScalarPolarization2D::TE_Ez,mat); }
    catch (const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected);
}

NOVA_TEST(helmholtz_2d_oblique_wave_on_perturbed_mesh) {
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,1,16,16);
    const C I(0,1);
    std::vector<bool> boundary(rect.mesh.nodes.size(),false);
    for (auto i:rect.boundary_nodes) boundary[i]=true;
    for (std::size_t i=0;i<rect.mesh.nodes.size();++i) if (!boundary[i]) {
        rect.mesh.nodes[i].x+=0.13/16*std::sin(static_cast<double>(i));
        rect.mesh.nodes[i].y+=0.13/16*std::cos(static_cast<double>(i));
    }
    auto exact=[&](std::size_t i) {
        const auto p=rect.mesh.nodes[i]; return std::exp(I*(1.2*p.x+1.6*p.y));
    };
    std::vector<HelmholtzDirichlet2D> bc;
    for (auto i:rect.boundary_nodes) bc.push_back({i,exact(i)});
    std::vector<HelmholtzMaterial2D> mat(rect.mesh.elements.size());
    const auto result=solve_helmholtz_2d(rect.mesh,2,ScalarPolarization2D::TE_Ez,mat,bc);
    for (std::size_t i=0;i<result.field.size();++i)
        NOVA_REQUIRE(std::abs(result.field[i]-exact(i))<0.006);
    NOVA_REQUIRE(result.relative_backward_error<1e-12);
}
NOVA_TEST(helmholtz_2d_dielectric_interface_te_tm) {
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,0.25,40,10);
    const C I(0,1);
    std::vector<HelmholtzMaterial2D> mat;
    for (auto e:rect.mesh.elements) {
        const double x=(rect.mesh.nodes[e[0]].x+rect.mesh.nodes[e[1]].x+
                        rect.mesh.nodes[e[2]].x)/3;
        mat.push_back({x<0.5?1.0:4.0,1.0});
    }
    for (auto pol:{ScalarPolarization2D::TE_Ez,ScalarPolarization2D::TM_Hz}) {
        const double q1=2, q2=4, a1=1, a2=pol==ScalarPolarization2D::TE_Ez?1:0.25;
        const double r=(a1*q1-a2*q2)/(a1*q1+a2*q2), t=1+r;
        NOVA_REQUIRE_NEAR(r*r+(a2*q2)/(a1*q1)*t*t,1,1e-14);
        std::vector<HelmholtzRobin2D> bc;
        for (auto e:rect.boundary_edges) {
            const auto x=rect.mesh.nodes[e[0]].x, y=rect.mesh.nodes[e[1]].x;
            if (x==0 && y==0) bc.push_back({e,-I*a1*q1,-2.0*I*a1*q1*std::exp(-I*q1*0.5)});
            if (x==1 && y==1) bc.push_back({e,-I*a2*q2,0});
        }
        const auto result=solve_helmholtz_2d(rect.mesh,2,pol,mat,{},bc);
        for (std::size_t i=0;i<result.field.size();++i) {
            const double x=rect.mesh.nodes[i].x-0.5;
            const C exact=x<=0?std::exp(I*q1*x)+r*std::exp(-I*q1*x):t*std::exp(I*q2*x);
            NOVA_REQUIRE(std::abs(result.field[i]-exact)<0.01);
        }
        NOVA_REQUIRE(result.relative_backward_error<1e-12);
    }
}
NOVA_TEST(helmholtz_2d_constant_volume_excitation) {
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,1,4,4);
    const C value(0.3,0.7), eps(2.0,0.5);
    std::vector<HelmholtzMaterial2D> mat(rect.mesh.elements.size(),{eps,1.0});
    std::vector<C> f(rect.mesh.elements.size(),-4.0*eps*value);
    const auto result=solve_helmholtz_2d(rect.mesh,2,ScalarPolarization2D::TE_Ez,mat,{},{},f);
    for (C u:result.field) NOVA_REQUIRE_NEAR(std::abs(u-value),0,1e-12);
    NOVA_REQUIRE_NEAR(std::accumulate(result.cell_loss.begin(),result.cell_loss.end(),0.0),
                      4.0*eps.imag()*std::norm(value),1e-12);
}

NOVA_TEST(helmholtz_2d_power_audit_dirichlet_reactions) {
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,1,1,1);
    const C value(0.3,0.7);
    std::vector<HelmholtzDirichlet2D> bc;
    for (auto i:rect.boundary_nodes) bc.push_back({i,value});
    for (auto pol:{ScalarPolarization2D::TE_Ez,ScalarPolarization2D::TM_Hz}) {
        std::vector<HelmholtzMaterial2D> mat(rect.mesh.elements.size(),
            {C(2,0.5),C(1,0.25)});
        const double expected=4.0*std::norm(value)*
            (pol==ScalarPolarization2D::TE_Ez?0.5:0.25);
        const auto result=solve_helmholtz_2d(rect.mesh,2,pol,mat,bc);
        NOVA_REQUIRE_NEAR(result.material_loss,expected,1e-12);
        NOVA_REQUIRE_NEAR(result.dirichlet_input,expected,1e-12);
        NOVA_REQUIRE_NEAR(result.excitation_input,0,1e-14);
        NOVA_REQUIRE_NEAR(result.robin_outflow,0,1e-14);
        NOVA_REQUIRE_NEAR(result.power_balance_residual,0,1e-12);
    }
}
NOVA_TEST(helmholtz_2d_power_audit_volume_excitation) {
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,1,4,4);
    const C value(0.3,0.7), eps(2.0,0.5);
    std::vector<HelmholtzMaterial2D> mat(rect.mesh.elements.size(),{eps,1.0});
    std::vector<C> f(rect.mesh.elements.size(),-4.0*eps*value);
    const auto result=solve_helmholtz_2d(rect.mesh,2,ScalarPolarization2D::TE_Ez,mat,{},{},f);
    const double expected=4.0*eps.imag()*std::norm(value);
    NOVA_REQUIRE_NEAR(result.excitation_input,expected,1e-12);
    NOVA_REQUIRE_NEAR(result.material_loss,expected,1e-12);
    NOVA_REQUIRE_NEAR(result.dirichlet_input,0,1e-14);
    NOVA_REQUIRE_NEAR(result.robin_outflow,0,1e-14);
    NOVA_REQUIRE_NEAR(result.power_balance_residual,0,1e-12);
}

NOVA_TEST(helmholtz_2d_sparse_dense_equivalence) {
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,0.5,8,4);
    const C value(0.3,0.7), I(0,1);
    std::vector<HelmholtzMaterial2D> mat(rect.mesh.elements.size(),
        {C(2,0.5),C(1,0.25)});
    std::vector<HelmholtzDirichlet2D> bc;
    std::vector<HelmholtzRobin2D> robin;
    for (auto i:rect.boundary_nodes) if (rect.mesh.nodes[i].x==0)
        bc.push_back({i,value});
    for (auto e:rect.boundary_edges)
        if (rect.mesh.nodes[e[0]].x==1 && rect.mesh.nodes[e[1]].x==1)
            robin.push_back({e,-I*2.0,C(0.2,0.1)});
    std::vector<C> sources(rect.mesh.elements.size(),C(0.4,-0.3));
    HelmholtzSolverOptions2D options;
    options.linear_solver=HelmholtzLinearSolver2D::SparseGmres;
    for (auto pol:{ScalarPolarization2D::TE_Ez,ScalarPolarization2D::TM_Hz}) {
        const auto dense=solve_helmholtz_2d(rect.mesh,2,pol,mat,bc,robin,sources);
        const auto sparse=solve_helmholtz_2d(rect.mesh,2,pol,mat,bc,robin,sources,options);
        for (std::size_t i=0;i<dense.field.size();++i)
            NOVA_REQUIRE(std::abs(dense.field[i]-sparse.field[i])<1e-9);
        NOVA_REQUIRE_NEAR(sparse.material_loss,dense.material_loss,1e-9);
        NOVA_REQUIRE_NEAR(sparse.power_balance_residual,0,1e-9);
        NOVA_REQUIRE(sparse.linear_iterations>0);
        NOVA_REQUIRE(sparse.matrix_nonzeros==dense.matrix_nonzeros);
    }
}
NOVA_TEST(helmholtz_2d_sparse_above_dense_node_limit) {
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,0.25,64,16);
    const C value(0.3,0.7), eps(2.0,0.5);
    std::vector<HelmholtzMaterial2D> mat(rect.mesh.elements.size(),{eps,1.0});
    std::vector<C> sources(rect.mesh.elements.size(),-4.0*eps*value);
    std::vector<HelmholtzDirichlet2D> bc;
    for (auto i:rect.boundary_nodes) bc.push_back({i,value});
    HelmholtzSolverOptions2D options;
    options.linear_solver=HelmholtzLinearSolver2D::SparseGmres;
    const auto result=solve_helmholtz_2d(rect.mesh,2,ScalarPolarization2D::TE_Ez,mat,bc,{},sources,options);
    NOVA_REQUIRE(result.field.size()>1024);
    NOVA_REQUIRE(result.matrix_nonzeros<8*result.field.size());
    for (C u:result.field) NOVA_REQUIRE(std::abs(u-value)<1e-8);
    NOVA_REQUIRE_NEAR(result.material_loss,eps.imag()*std::norm(value),1e-8);
    NOVA_REQUIRE_NEAR(result.power_balance_residual,0,1e-8);
    NOVA_REQUIRE(result.relative_backward_error<1e-10);
}
