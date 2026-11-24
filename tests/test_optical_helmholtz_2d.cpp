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
    auto rect=nova::fem::make_rectangle_mesh_2d(0,1,0,0.25,nx,2);
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
