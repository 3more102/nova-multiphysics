#include "nova/physics/oblique_optics.hpp"
#include <array>
#include <iostream>
#include <stdexcept>
using namespace nova::optics;
void near(double x,double y,double tol) {
    if (!std::isfinite(x) || std::abs(x-y)>tol) throw std::runtime_error("reference mismatch");
}
int main() {
    constexpr double l=1.55e-6;
    for (auto p:{Polarization::TE,Polarization::TM}) {
        auto n=oblique_incidence_stack(l,1,1.5,0,p,{});
        near(n.reflectance,0.04,1e-12); near(n.transmittance,0.96,1e-12);
        auto tir=oblique_incidence_stack(l,1.5,1,1.0,p,{});
        near(tir.reflectance,1,1e-12); near(tir.transmittance,0,1e-12);
        const std::array<Layer,1> gap{{{1.0,l/10}}};
        auto ft=oblique_incidence_stack(l,1.5,1.5,1.0,p,gap);
        if (!(ft.transmittance>0 && ft.transmittance<1)) return 1;
        near(ft.absorptance,0,1e-12);
        const std::array<Layer,1> film{{{Complex{1.5,0.5},1e-3}}};
        auto loss=oblique_incidence_stack(l,1,1,0.4,p,film);
        near(loss.transmittance,0,1e-12);
        if (!(loss.absorptance>0)) return 1;
        const std::array<Layer,1> match{{{std::sqrt(1.5),l/(4*std::sqrt(1.5))}}};
        near(oblique_incidence_stack(l,1,1.5,0,p,match).reflectance,0,1e-12);
        bool rejected=false;
        try { (void)oblique_incidence_stack(l,1,1,-0.1,p,{}); }
        catch (const std::invalid_argument&) { rejected=true; }
        if (!rejected) return 1;
    }
    auto brewster=oblique_incidence_stack(l,1,1.5,std::atan(1.5),Polarization::TM,{});
    near(brewster.reflectance,0,1e-12); near(brewster.transmittance,1,1e-12);
    const double a=0.6, b=std::asin(std::sin(a)/1.5);
    auto te=oblique_incidence_stack(l,1,1.5,a,Polarization::TE,{});
    near(te.reflectance,std::pow((std::cos(a)-1.5*std::cos(b))/(std::cos(a)+1.5*std::cos(b)),2),1e-12);
    auto tm=oblique_incidence_stack(l,1,1.5,a,Polarization::TM,{});
    near(tm.reflectance,std::pow((1.5*std::cos(a)-std::cos(b))/(1.5*std::cos(a)+std::cos(b)),2),1e-12);
    std::cout << "Oblique optics reference cases passed\n";
}
