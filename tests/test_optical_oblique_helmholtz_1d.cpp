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
// Independent exact critical-angle solution in a uniform exit-matched slab.
// The scalar field is identically 2 and the outgoing flux vanishes.
void critical_exit_limit() {
    const double ni=1.5, theta=0.7, nt=ni*std::sin(theta);
    std::vector<double> nodes;
    for(std::size_t j=0;j<=80;++j)
        nodes.push_back(0.37*wavelength*static_cast<double>(j)/80.0);
    for(auto pol:{Polarization::te,Polarization::tm}) {
        const bool tm=pol==Polarization::tm;
        for(auto material:{Complex{nt*nt,0},Complex{2.1,0.3}}) {
            const std::vector<Complex> eps(80,material);
            const auto exact=solve_oblique_helmholtz_1d(wavelength,ni,nt,
                theta,pol,nodes,eps);
            require(exact.transmittance==0,"critical exit has zero normal flux");
            require(std::abs(exact.power_balance_residual)<1e-9,
                "critical exit power balance");
            if(material.imag()==0) {
                require(std::abs(exact.reflection-(tm?-1.0:1.0))<1e-9,
                    "critical exit exact reflected amplitude");
                require(std::abs(exact.transmission-(tm?0.0:2.0))<1e-9,
                    "critical exit tangential electric amplitude");
                require(exact.absorptance==0,"lossless critical exit absorption");
                for(auto u:exact.scalar_field)
                    require(std::abs(u-2.0)<1e-9,"critical exit constant scalar field");
            } else {
                require(exact.absorptance>0 && exact.absorptance<1,
                    "lossy critical exit absorbs power");
                for(double a:exact.cell_absorptance)
                    require(a>=0,"critical exit passive cell loss");
            }
            for(double side:{-1.0,1.0}) {
                const auto near=solve_oblique_helmholtz_1d(wavelength,ni,nt,
                    theta+side*1e-9,pol,nodes,eps);
                require(std::abs(near.reflection-exact.reflection)<5e-4 &&
                    std::abs(near.transmission-exact.transmission)<5e-4 &&
                    std::abs(near.absorptance-exact.absorptance)<5e-4,
                    "continuous amplitude and absorption across exit critical angle");
                require(std::abs(near.power_balance_residual)<1e-9,
                    "near-critical exit power balance");
                require(side<0?near.transmittance>0:near.transmittance==0,
                    "propagating versus evanescent exit flux");
            }
        }
    }
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
        critical_exit_limit();
        std::cout<<"NOVA polarized oblique Helmholtz FEM tests passed\n";
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
