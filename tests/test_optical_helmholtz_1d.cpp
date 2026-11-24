#include "nova/physics/optical_helmholtz_1d.hpp"
#include "nova/physics/optics.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
using namespace nova::optics;
void require(bool b, const char* msg) { if (!b) throw std::runtime_error(msg); }
constexpr double wavelength=1.55e-6;
Helmholtz1DResult slab(std::size_t count, Complex index, double length) {
    std::vector<double> nodes(count+1);
    for (std::size_t j=0; j<=count; ++j)
        nodes[j]=length*static_cast<double>(j)/static_cast<double>(count);
    const std::vector<Complex> eps(count,index*index);
    return solve_helmholtz_1d(wavelength,1.0,1.0,nodes,eps);
}
}
int main() {
    try {
        const Complex exact = std::exp(Complex{0,2*std::numbers::pi_v<double>*0.73});
        const auto coarse=slab(40,1.0,0.73*wavelength);
        const auto fine=slab(80,1.0,0.73*wavelength);
        require(std::abs(fine.transmission-exact)<0.30*std::abs(coarse.transmission-exact),
                "second-order uniform wave convergence");
        require(fine.reflectance<1e-6 && std::abs(fine.power_balance_residual)<1e-11,
                "lossless uniform power");
        for (Complex index : {Complex{1.7,0},Complex{1.7,0.25},Complex{0.2,3.0}}) {
            const double d=0.37*wavelength;
            const std::array<Layer,1> layers{{{index,d}}};
            const auto reference=normal_incidence_stack(wavelength,1.0,1.0,layers);
            const auto fem=slab(400,index,d);
            require(std::abs(fem.reflection-reference.reflection)<2e-4,"slab complex reflection");
            require(std::abs(fem.transmission-reference.transmission)<2e-4,"slab complex transmission");
            require(std::abs(fem.absorptance-reference.absorptance)<2e-4,"slab absorption");
            require(std::abs(fem.power_balance_residual)<1e-10,"discrete energy balance");
            for (double a : fem.cell_absorptance) require(a>=0,"passive cell loss");
        }
        // A deliberately coarse element forces an adjacent pivot: diagonal real
        // part vanishes at h=sqrt(3), while the off-diagonal exceeds the end load.
        const std::array<double,2> x{0,wavelength*std::sqrt(3.0)/(2*std::numbers::pi_v<double>)};
        const std::array<Complex,1> e{1.0};
        const auto pivot=solve_helmholtz_1d(wavelength,0.1,0.1,x,e);
        require(std::abs(pivot.power_balance_residual)<1e-12,"pivoted solve power");
        bool rejected=false;
        const std::array<Complex,1> gain{Complex{1,-0.1}};
        try { (void)solve_helmholtz_1d(wavelength,1,1,x,gain); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"reject gain");
        rejected=false;
        const std::array<double,2> repeated{0,0};
        try { (void)solve_helmholtz_1d(wavelength,1,1,repeated,e); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"reject zero-length element");
        std::cout<<"NOVA optical Helmholtz FEM reference tests passed\n";
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
