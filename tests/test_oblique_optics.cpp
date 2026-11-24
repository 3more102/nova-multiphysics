#include "nova/physics/oblique_optics.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void near(double a,double b,double tol,const char* msg) {
    if (!std::isfinite(a) || std::abs(a-b)>tol) throw std::runtime_error(msg);
}
void require(bool ok,const char* msg) {
    if (!ok) throw std::runtime_error(msg);
}
}

int main() {
    using namespace nova::optics;
    constexpr double l=1.55e-6;
    constexpr double k=2*std::numbers::pi_v<double>/l;
    try {
        const std::array<Layer,2> stack{{{{2.0,0.05},l/8},{{1.4,0},l/3}}};
        const auto normal=normal_incidence_stack(l,1,1.5,stack);
        for (auto pol : {Polarization::TE,Polarization::TM}) {
            const auto a=oblique_incidence_stack(l,1,1.5,0,pol,stack);
            near(std::abs(a.reflection-normal.reflection),0,1e-12,"normal r agreement");
            near(std::abs(a.transmission-normal.transmission),0,1e-12,"normal t agreement");
            near(a.absorptance,normal.absorptance,1e-12,"normal A agreement");
        }
        const double theta=std::atan(1.5);
        const auto brew=oblique_incidence_stack(l,1,1.5,theta,Polarization::TM,{});
        near(brew.reflectance,0,1e-13,"Brewster TM");
        near(brew.transmittance,1,1e-12,"Brewster T");
        const double qi=std::cos(theta), qe=std::sqrt(2.25-std::pow(std::sin(theta),2));
        const auto te=oblique_incidence_stack(l,1,1.5,theta,Polarization::TE,{});
        near(te.reflectance,std::pow((qi-qe)/(qi+qe),2),1e-13,"oblique TE Fresnel");
        near(te.reflectance+te.transmittance,1,1e-12,"TE power");
        const double tir=std::numbers::pi_v<double>/3;
        for (auto pol : {Polarization::TE,Polarization::TM}) {
            const auto a=oblique_incidence_stack(l,1.5,1,tir,pol,{});
            near(a.reflectance,1,1e-12,"TIR R");
            near(a.transmittance,0,1e-14,"TIR normal power");
            near(a.absorptance,0,1e-12,"TIR absorption");
        }
        const double s=1.5*std::sin(tir);
        const double q=1.5*std::cos(tir), alpha=std::sqrt(s*s-1);
        const double d=l/5;
        const std::array<Layer,1> gap{{{1.0,d}}};
        const double exact=1/(1+std::pow((q*q+alpha*alpha)/(2*q*alpha),2)
            *std::pow(std::sinh(k*alpha*d),2));
        const auto tunnel=oblique_incidence_stack(l,1.5,1.5,tir,Polarization::TE,gap);
        near(tunnel.transmittance,exact,1e-12,"frustrated TIR analytical T");
        near(tunnel.reflectance+tunnel.transmittance,1,1e-12,"tunneling power");
        const std::array<Layer,1> metal{{{{0.2,4.0},1e-3}}};
        for(auto pol : {Polarization::TE,Polarization::TM}) {
            const auto a=oblique_incidence_stack(l,1,1,0.4,pol,metal);
            near(a.transmittance,0,1e-14,"thick metal T");
            require(a.absorptance>0 && a.reflectance<1,"thick metal absorption");
        }
        const std::array<Layer,1> zero{{{2.0,0.0}}};
        const auto z=oblique_incidence_stack(l,1,1.5,theta,Polarization::TE,zero);
        near(z.reflectance,te.reflectance,1e-13,"zero layer invariance");
        // Reciprocity at the same conserved transverse wavenumber.
        const double a0=0.5, a1=std::asin(std::sin(a0)/1.5);
        const std::array<Layer,2> reverse{{stack[1],stack[0]}};
        for(auto pol : {Polarization::TE,Polarization::TM}) {
            const auto f=oblique_incidence_stack(l,1,1.5,a0,pol,stack);
            const auto b=oblique_incidence_stack(l,1.5,1,a1,pol,reverse);
            near(f.transmittance,b.transmittance,1e-12,"reciprocity T");
        }
        bool rejected=false;
        try { (void)oblique_incidence_stack(l,1,1,-0.1,Polarization::TE,{}); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"negative angle rejection");
        rejected=false;
        const std::array<Layer,1> gain{{{{1,-0.1},l}}};
        try { (void)oblique_incidence_stack(l,1,1,0.2,Polarization::TM,gain); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"gain rejection");
        std::cout<<"NOVA oblique TE/TM reference cases passed\n";
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
