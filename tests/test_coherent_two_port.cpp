#include "nova/physics/coherent_two_port.hpp"
#include <array>
#include <iostream>
#include <stdexcept>

namespace {
void near(double a,double b,double tol,const char* msg) {
    if (!std::isfinite(a) || std::abs(a-b)>tol) throw std::runtime_error(msg);
}
}
int main() {
    using namespace nova::optics;
    constexpr double lambda=1.55e-6;
    try {
        for (auto pol : {Polarization::te,Polarization::tm}) {
            const auto bare=coherent_two_port(lambda,1.0,1.7,0.4,pol,{});
            near(bare.minimum_absorptance,0,1e-12,"lossless minimum");
            near(bare.maximum_absorptance,0,1e-12,"lossless maximum");
            near(std::abs(bare.s12-bare.s21),0,1e-12,"reciprocity");
            near(coherent_power(bare,{1,0},{0.3,0.7}).absorbed_power,0,1e-12,"unitarity");
            const std::array<Layer,3> asym{{{{2.0,0.1},lambda/5},
                {{1.4,0.0},lambda/7},{{1.8,0.3},lambda/3}}};
            const auto s=coherent_two_port(lambda,1.0,1.5,0.5,pol,asym);
            const auto one=oblique_incidence_stack(lambda,1,1.5,0.5,pol,asym);
            near(coherent_power(s,1,0).absorbed_power,one.absorptance,1e-12,"one-port parity");
            for(int j=0;j<64;++j) {
                const Complex phase=std::polar(1.0,2*std::numbers::pi_v<double>*j/64);
                const double a=coherent_power(s,1,phase).absorbed_power/2;
                if(a<s.minimum_absorptance-1e-12 || a>s.maximum_absorptance+1e-12)
                    throw std::runtime_error("coherent bounds");
            }
        }
        const Complex n{2.0,0.3};
        const double d=lambda/5;
        const std::array<Layer,1> film{{{n,d}}};
        const auto slab=coherent_two_port(lambda,1,1,0,Polarization::te,film);
        const Complex rho=(1.0-n)/(1.0+n);
        const Complex p=std::exp(Complex{0,2*std::numbers::pi_v<double>*d/lambda}*n);
        // Independent even/odd Fabry-Perot scattering eigenvalues.
        const double even=1-std::norm((rho+p)/(1.0+rho*p));
        const double odd=1-std::norm((rho-p)/(1.0-rho*p));
        near(coherent_power(slab,1,1).absorbed_power/2,even,1e-12,"even eigenchannel");
        near(coherent_power(slab,1,-1).absorbed_power/2,odd,1e-12,"odd eigenchannel");
        near(slab.minimum_absorptance,std::min(even,odd),1e-12,"minimum eigenvalue");
        near(slab.maximum_absorptance,std::max(even,odd),1e-12,"maximum eigenvalue");
        if(std::abs(even-odd)<0.01) throw std::runtime_error("phase contrast fixture");
        near(coherent_power(slab,0,0).absorbed_power,0,0,"zero excitation");
        bool rejected=false;
        try { (void)coherent_two_port(lambda,1.5,1,1.0,Polarization::te,{}); }
        catch(const std::domain_error&) { rejected=true; }
        if(!rejected) throw std::runtime_error("evanescent port rejection");
        std::cout<<"Coherent optical two-port tests passed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
