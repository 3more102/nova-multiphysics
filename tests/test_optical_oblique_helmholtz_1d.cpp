#include "nova/physics/optical_helmholtz_1d.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
using namespace nova::optics;
constexpr double wavelength=1.55e-6;
void require(bool ok,const char* msg) {
    if(!ok) throw std::runtime_error(msg);
}
ObliqueHelmholtz1DResult run(double ni,double nt,double angle,Polarization pol,
    std::span<const Layer> layers,std::size_t m) {
    std::vector<double> nodes{0};
    std::vector<Complex> eps;
    for(const auto& l:layers) {
        const double start=nodes.back();
        for(std::size_t j=1;j<=m;++j) {
            nodes.push_back(start+l.thickness_m*static_cast<double>(j)/static_cast<double>(m));
            eps.push_back(l.index*l.index);
        }
    }
    return solve_oblique_helmholtz_1d(wavelength,ni,nt,angle,pol,nodes,eps);
}
void reciprocity_and_normal_limit(std::span<const Layer> layers) {
    constexpr double ni=1.3, nt=1.8, angle=0.62;
    const double reverse_angle=std::asin(ni*std::sin(angle)/nt);
    const std::vector<Layer> reversed(layers.rbegin(),layers.rend());
    for(auto pol:{Polarization::te,Polarization::tm}) {
        const auto forward=run(ni,nt,angle,pol,layers,200);
        const auto reverse=run(nt,ni,reverse_angle,pol,reversed,200);
        const double qi=ni*std::cos(angle), qt=nt*std::cos(reverse_angle);
        const double yi=pol==Polarization::te?qi:ni*ni/qi;
        const double yt=pol==Polarization::te?qt:nt*nt/qt;
        require(std::abs(forward.transmittance-reverse.transmittance)<1e-9,
            "reciprocal transmitted power on reversed mesh");
        require(std::abs(std::sqrt(yt/yi)*forward.transmission-
            std::sqrt(yi/yt)*reverse.transmission)<1e-9,
            "reciprocal complex power-normalized transmission");
        require(std::abs(reverse.power_balance_residual)<1e-9,
            "reverse incidence energy balance");
    }
    std::vector<double> nodes{0};
    std::vector<Complex> eps;
    for(const auto& layer:layers) {
        const double start=nodes.back();
        for(std::size_t j=1;j<=80;++j) {
            nodes.push_back(start+layer.thickness_m*static_cast<double>(j)/80.0);
            eps.push_back(layer.index*layer.index);
        }
    }
    const auto normal=solve_helmholtz_1d(wavelength,ni,nt,nodes,eps);
    const auto te=solve_oblique_helmholtz_1d(wavelength,ni,nt,0,
        Polarization::te,nodes,eps);
    require(std::abs(normal.reflection-te.reflection)<1e-12,
        "normal and TE zero-angle reflection agree");
    require(std::abs(normal.transmission-te.transmission)<1e-12,
        "normal and TE zero-angle transmission agree");
    require(normal.electric_field.size()==te.scalar_field.size(),
        "normal and TE field sizes agree");
    require(normal.cell_absorptance.size()==te.cell_absorptance.size(),
        "normal and TE loss sizes agree");
    for(std::size_t j=0;j<nodes.size();++j)
        require(std::abs(normal.electric_field[j]-te.scalar_field[j])<1e-12,
            "normal and TE nodal fields agree");
    for(std::size_t j=0;j<eps.size();++j)
        require(std::abs(normal.cell_absorptance[j]-te.cell_absorptance[j])<1e-12,
            "normal and TE cell loss agree");
}
// Independent lossless symmetric barrier benchmark. Inside the barrier
// the longitudinal wavenumber is imaginary, but both ports propagate.
void frustrated_total_internal_reflection() {
    constexpr double n=1.5, barrier=1.0, theta=0.95;
    const double b=n*std::sin(theta);
    const double q=std::sqrt(n*n-b*b);
    const double kappa=std::sqrt(b*b-barrier*barrier);
    for(auto pol:{Polarization::te,Polarization::tm}) {
        const double g=pol==Polarization::te?q:q/(n*n);
        const double v=pol==Polarization::te?kappa:kappa/(barrier*barrier);
        double previous=1.0;
        for(double d:{0.05,0.25,0.75,1.5}) {
            const double s=std::sinh(2*std::numbers::pi*kappa*d);
            const double factor=(g*g+v*v)/(2*g*v);
            const double exact=1/(1+factor*factor*s*s);
            const std::array<Layer,1> layers{{{barrier,d*wavelength}}};
            const auto coarse=run(n,n,theta,pol,layers,100);
            const auto fine=run(n,n,theta,pol,layers,400);
            require(fine.transmittance>0 && fine.transmittance<previous,
                "finite evanescent barrier transmits with thickness decay");
            require(std::abs(fine.transmittance-exact)<3e-7,
                "independent frustrated TIR transmission");
            require(std::abs(fine.reflectance-(1-exact))<3e-7,
                "independent frustrated TIR reflection");
            require(std::abs(fine.transmittance-exact)<0.09*
                std::abs(coarse.transmittance-exact)+1e-10,
                "frustrated TIR mesh convergence");
            require(std::abs(fine.absorptance)<1e-12,
                "lossless evanescent barrier does not absorb");
            require(std::abs(fine.power_balance_residual)<1e-9,
                "frustrated TIR power conservation");
            previous=fine.transmittance;
        }
    }
}
void compare(double ni,double nt,double angle,std::span<const Layer> layers) {
    for(auto pol:{Polarization::te,Polarization::tm}) {
        const auto exact=oblique_incidence_stack(wavelength,ni,nt,angle,pol,layers);
        const auto coarse=run(ni,nt,angle,pol,layers,200);
        const auto fine=run(ni,nt,angle,pol,layers,400);
        const double err=std::abs(fine.transmission-exact.transmission);
        require(err<1e-4,"oblique complex transmission");
        require(std::abs(fine.reflection-exact.reflection)<1e-4,"oblique complex reflection");
        require(std::abs(fine.reflectance-exact.reflectance)<1e-4,"oblique R");
        require(std::abs(fine.transmittance-exact.transmittance)<1e-4,"oblique T");
        require(std::abs(fine.absorptance-exact.absorptance)<1e-4,"oblique A");
        require(std::abs(fine.power_balance_residual)<1e-9,"oblique energy balance");
        require(err<0.30*std::abs(coarse.transmission-exact.transmission)+1e-10,
            "oblique second-order transmission convergence");
        for(double a:fine.cell_absorptance) require(a>=0,"passive cell loss");
        if(ni*std::sin(angle)>nt) require(fine.transmittance==0,"evanescent exit power");
    }
}
}
int main() {
    try {
        const std::array<Layer,3> lossy{{
            {1.7,0.18*wavelength},{{0.25,2.0},0.07*wavelength},
            {{2.2,0.08},0.31*wavelength}}};
        compare(1.3,1.8,0.62,lossy);
        reciprocity_and_normal_limit(lossy);
        const std::array<Layer,1> tir{{{1.2,0.25*wavelength}}};
        compare(1.5,1.0,0.95,tir);
        frustrated_total_internal_reflection();
        const std::array<Layer,1> brewster{{{1.5,0.2*wavelength}}};
        const double angle=std::atan(1.5);
        compare(1.0,1.5,angle,brewster);
        require(run(1,1.5,angle,Polarization::tm,brewster,400).reflectance<1e-9,
            "TM Brewster reflection vanishes");
        const std::array<Layer,1> normal{{{{1.7,0.25},0.37*wavelength}}};
        compare(1,1,0,normal);
        const std::array<double,2> x{0,wavelength/10};
        const std::array<Complex,1> e{{{0,0}}};
        bool rejected=false;
        try { (void)solve_oblique_helmholtz_1d(wavelength,1,1,0.2,Polarization::tm,x,e); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"reject zero TM permittivity");
        rejected=false;
        try { (void)solve_oblique_helmholtz_1d(wavelength,1,1,-0.2,Polarization::te,x,e); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"reject negative incidence angle");
        std::cout<<"NOVA polarized oblique Helmholtz FEM tests passed\n";
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
