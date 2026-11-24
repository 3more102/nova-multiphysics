#include "nova/physics/optics_bragg.hpp"
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {
void near(double a, double b, double tol) {
    if (!std::isfinite(a) || std::abs(a-b)>tol)
        throw std::runtime_error("Bragg reference mismatch");
}
}
int main() {
    using nova::optics::uniform_bragg_grating;
    constexpr double tol=1e-12;
    try {
        const auto empty=uniform_bragg_grating(3,7,0);
        near(empty.reflectance,0,tol); near(empty.transmittance,1,tol);
        const auto uncoupled=uniform_bragg_grating(3,0,2);
        near(std::abs(uncoupled.transmission-std::polar(1.0,6.0)),0,tol);
        const auto resonance=uniform_bragg_grating(0,2,1);
        near(resonance.reflectance,std::pow(std::tanh(2.0),2),tol);
        near(resonance.transmittance,1/std::pow(std::cosh(2.0),2),tol);
        near(resonance.reflection.imag(),std::tanh(2.0),tol);
        const auto edge=uniform_bragg_grating(2,2,1);
        near(edge.reflectance,0.8,tol); near(edge.transmittance,0.2,tol);
        for(double d : {2.0-1e-10,2.0+1e-10})
            near(uniform_bragg_grating(d,2,1).reflectance,0.8,1e-9);
        const auto null=uniform_bragg_grating(std::hypot(2.0,std::numbers::pi),2,1);
        near(null.reflectance,0,tol); near(null.transmittance,1,tol);
        const auto strong=uniform_bragg_grating(0,1000,1);
        near(strong.reflectance,1,tol); near(strong.transmittance,0,tol);
        for(int i=-200;i<=200;++i) {
            const double d=0.05*i;
            const auto p=uniform_bragg_grating(d,2,1);
            const auto q=uniform_bragg_grating(-d,-2,1);
            near(p.reflectance+p.transmittance,1,tol);
            near(std::abs(p.reflection-std::conj(q.reflection)),0,tol);
            near(std::abs(p.transmission-std::conj(q.transmission)),0,tol);
        }
        bool rejected=false;
        try { (void)uniform_bragg_grating(1,1,-1); }
        catch(const std::invalid_argument&) { rejected=true; }
        if(!rejected) throw std::runtime_error("negative length accepted");
        std::cout << "NOVA uniform Bragg reference passed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
