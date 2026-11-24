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
        const std::array<Layer,1> tir{{{1.2,0.25*wavelength}}};
        compare(1.5,1.0,0.95,tir);
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
