#include "nova/physics/kerr_cavity.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool yes, const char* message) {
    if (!yes) throw std::runtime_error(message);
}
void near(double a, double b, double tol=1e-12) {
    require(std::isfinite(a) && std::abs(a-b)<=tol, "reference mismatch");
}
}
int main() {
    using namespace nova::optics;
    using S=KerrCavityStability;
    try {
        const auto zero=kerr_cavity_steady_states(3,0);
        require(zero.size()==1 && zero[0].intensity==0 && zero[0].stability==S::stable,"zero drive");
        const auto on=kerr_cavity_steady_states(0,2);
        require(on.size()==1 && on[0].stability==S::stable,"monostable resonance");
        near(on[0].intensity,1);
        const auto three=kerr_cavity_steady_states(3,3);
        require(three.size()==3,"three distinct roots");
        near(three[0].intensity,(3-std::sqrt(5.0))/2);
        near(three[1].intensity,(3+std::sqrt(5.0))/2);
        near(three[2].intensity,3);
        require(three[0].stability==S::stable && three[1].stability==S::unstable && three[2].stability==S::stable,"branch stability");
        const auto fold=kerr_cavity_steady_states(2,2);
        require(fold.size()==2 && fold[0].stability==S::marginal && fold[1].stability==S::stable,"exact fold");
        near(fold[0].intensity,1); near(fold[1].intensity,2);
        require(kerr_cavity_steady_states(2,2-1e-8).size()==3,"below fold");
        require(kerr_cavity_steady_states(2,2+1e-8).size()==1,"above fold");
        for (double d : {-1000.,-3.,0.,1.,2.,3.,10.,1000.}) {
            for (double p : {1e-12,0.01,1.,2.,3.,10.,100.,1e6}) {
                const auto states=kerr_cavity_steady_states(d,p);
                require(!states.empty() && states.size()<=3,"root count");
                double prev=-1;
                for (const auto& state : states) {
                    const double x=state.intensity;
                    require(x>prev && x<=p,"root bounds/order");
                    near(x*(1+(d-x)*(d-x))/p,1,2e-10);
                    prev=x;
                }
            }
        }
        for (double p : {-1.,1e-13,1e7,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
            bool rejected=false;
            try { (void)kerr_cavity_steady_states(0,p); } catch (const std::invalid_argument&) { rejected=true; }
            require(rejected,"reject invalid drive");
        }
        std::cout << "Kerr cavity analytical benchmarks and 64 sweep cases passed\n";
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
