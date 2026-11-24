#include "nova/physics/polarization.hpp"
#include <iostream>
#include <limits>

int main() {
    using namespace nova::optics;
    using C = std::complex<double>;
    constexpr double lambda=1.55e-6, pi=std::numbers::pi_v<double>;
    const JonesField x{C{1},C{0}};
    const auto power=[](const JonesField& e){return std::norm(e[0])+std::norm(e[1]);};
    const auto near=[](C a,C b){if (!(std::abs(a-b)<1e-12)) throw std::runtime_error("Jones reference mismatch");};
    try {
        near(apply_jones(jones_stack(lambda,{}),x)[0],1);
        const std::array<JonesLayer,1> half{{{1.0,1.5,lambda,pi/4}}};
        const auto h=apply_jones(jones_stack(lambda,half),x);
        near(h[0],0); near(std::norm(h[1]),1);
        const std::array<JonesLayer,1> quarter{{{1.0,1.5,lambda/2,pi/4}}};
        const auto q=apply_jones(jones_stack(lambda,quarter),x);
        near(std::norm(q[0]),0.5); near(q[1]/q[0],C{0,-1});
        const std::array<JonesLayer,1> loss{{{C{1,0.1},C{1,0.2},lambda,0}}};
        near(power(apply_jones(jones_stack(lambda,loss),x)),std::exp(-4*pi*0.1));
        const std::array<JonesLayer,2> seq{{half[0],loss[0]}}, rev{{loss[0],half[0]}};
        near(power(apply_jones(jones_stack(lambda,seq),x)),std::exp(-4*pi*0.2));
        near(power(apply_jones(jones_stack(lambda,rev),x)),std::exp(-4*pi*0.1));
        const std::array<JonesLayer,1> iso{{{1.5,1.5,lambda,0.37}}};
        const auto i=apply_jones(jones_stack(lambda,iso),x);
        near(i[0],-1); near(i[1],0);
        auto bad=loss; bad[0].fast_index=C{1,-0.1};
        bool rejected=false;
        try {(void)jones_stack(lambda,bad);} catch (const std::invalid_argument&) {rejected=true;}
        if (!rejected) throw std::runtime_error("gain accepted");
        rejected=false;
        try {(void)jones_stack(std::numeric_limits<double>::quiet_NaN(),half);} catch (const std::invalid_argument&) {rejected=true;}
        if (!rejected) throw std::runtime_error("NaN accepted");
        std::cout<<"Jones polarization references passed\n";
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n'; return 1;}
}
