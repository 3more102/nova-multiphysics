#include "nova/physics/optics.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double a, double b, double tol, const char* message) {
    require(std::isfinite(a) && std::abs(a-b) <= tol, message);
}
}
int main() {
    using namespace nova::optics;
    constexpr double lambda = 1.55e-6;
    constexpr double k = 2.0 * std::numbers::pi_v<double> / lambda;
    try {
        const auto f = normal_incidence_stack(lambda, 1.0, 1.5, {});
        near(f.reflectance, 0.04, 1e-13, "Fresnel R");
        near(f.transmittance, 0.96, 1e-13, "Fresnel T");
        const std::array<Layer,1> match{{{std::sqrt(1.5), lambda/(4*std::sqrt(1.5))}}};
        const auto ar = normal_incidence_stack(lambda,1,1.5,match);
        near(ar.reflectance,0,1e-13,"quarter wave antireflection");
        near(ar.transmittance,1,1e-13,"AR T");
        std::vector<Layer> dbr;
        for (int i=0;i<8;++i) {
            dbr.push_back({2.0,lambda/8.0}); dbr.push_back({1.5,lambda/6.0});
        }
        const double q = std::pow(2.0/1.5,16);
        const auto br = normal_incidence_stack(lambda,1,1,dbr);
        near(br.reflectance,std::pow((q-1)/(q+1),2),1e-12,"DBR analytical R");
        near(br.reflectance+br.transmittance,1,1e-12,"DBR energy");
        const std::array<Layer,1> lossy{{{Complex{1.5,0.5}, 1e-3}}};
        const auto loss = normal_incidence_stack(lambda,1,1,lossy);
        near(loss.reflectance,std::norm((1.0-lossy[0].index)/(1.0+lossy[0].index)),1e-12,"thick absorbing layer");
        near(loss.transmittance,0,1e-12,"thick layer T");
        require(loss.absorptance>0,"passive absorption");
        const auto eps = drude_permittivity(2.0,1.0,6.0,0.0);
        near(eps.real(),-8.0,1e-14,"lossless Drude");
        require(drude_permittivity(2,1,6,0.1).imag()>0,"passive Drude");
        near(passive_index({-4,0}).imag(),2,1e-14,"passive square root");
        const auto spp = planar_surface_plasmon(lambda,1,{-10,0});
        near(spp.beta_per_m.real()/k,std::sqrt(10.0/9.0),1e-13,"SPP dispersion");
        require(std::isinf(spp.intensity_propagation_length_m),"lossless SPP");
        const auto sp = planar_surface_plasmon(lambda,1,{-10,1});
        require(sp.beta_per_m.imag()>0 && sp.intensity_propagation_length_m>0,"lossy SPP");
        require(std::abs(sp.dielectric_decay_per_m + sp.metal_decay_per_m/Complex{-10,1})/k < 1e-12,"TM interface condition");
        const auto kr = kerr_phase(lambda,2e-19,1e12,0.01);
        near(kr.index_shift,2e-7,1e-19,"Kerr index");
        near(kr.phase_shift_rad,k*2e-9,1e-13,"Kerr phase");
        near(kerr_phase(lambda,-2e-19,1e12,0.01).index_shift,-2e-7,1e-19,"negative Kerr");
        bool rejected=false;
        try { (void)planar_surface_plasmon(lambda,1,{-0.5,0}); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"reject unsupported SPP");
        rejected=false;
        try { (void)passive_index({1,-1}); } catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"reject gain");
        std::cout << "NOVA optics analytical reference cases passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
