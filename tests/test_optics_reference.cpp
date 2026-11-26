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
        // Zero-width layers must be exactly equivalent to their removal, even
        // when their interface reflection rounds to +/-1 at extreme contrast.
        const std::array<Layer,2> empty_stack{{{1e30,0.0},{1e-30,0.0}}};
        const auto empty_result = normal_incidence_stack(lambda,1.0,1.5,empty_stack);
        require(empty_result.reflection == f.reflection &&
                empty_result.transmission == f.transmission,"zero-width interface identity");
        const std::array<Layer,1> finite_film{{{Complex{2.0,0.1},lambda/7}}};
        const std::array<Layer,3> padded_film{{{1e30,0.0},finite_film[0],{1e-30,0.0}}};
        const auto unpadded = normal_incidence_stack(lambda,1.0,1.5,finite_film);
        const auto padded = normal_incidence_stack(lambda,1.0,1.5,padded_film);
        require(padded.reflection == unpadded.reflection &&
                padded.transmission == unpadded.transmission &&
                padded.absorptance == unpadded.absorptance,"zero-width lossy-film identity");
        bool invalid_zero_layer = false;
        try {
            const std::array<Layer,1> bad{{{Complex{1.0,-0.1},0.0}}};
            (void)normal_incidence_stack(lambda,1.0,1.5,bad);
        } catch (const std::invalid_argument&) { invalid_zero_layer = true; }
        require(invalid_zero_layer,"zero-width layers retain material validation");
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

        // Oblique incidence: independent single-interface Fresnel and power checks.
        const double theta = std::numbers::pi_v<double>/6.0;
        const double ci = std::cos(theta), ct = std::sqrt(1-std::pow(std::sin(theta)/1.5,2));
        const auto te = oblique_incidence_stack(lambda,1,1.5,theta,Polarization::te,{});
        const auto tm = oblique_incidence_stack(lambda,1,1.5,theta,Polarization::tm,{});
        near(te.reflectance,std::pow((ci-1.5*ct)/(ci+1.5*ct),2),1e-13,"oblique TE Fresnel");
        near(tm.reflectance,std::pow((1.5*ci-ct)/(1.5*ci+ct),2),1e-13,"oblique TM Fresnel");
        const auto brewster = oblique_incidence_stack(lambda,1,1.5,std::atan(1.5),Polarization::tm,{});
        near(brewster.reflectance,0,1e-13,"Brewster reflection");
        near(brewster.transmittance,1,1e-13,"Brewster transmission");
        const std::array<Layer,1> gap{{{1.0,lambda/10}}};
        const std::array<Layer,1> wide_gap{{{1.0,lambda*10}}};
        const std::array<Layer,1> film{{{2.0,lambda/7}}};
        const std::array<Layer,1> zero{{{1.0,0}}};
        for (const auto pol : {Polarization::te,Polarization::tm}) {
            const auto parity = oblique_incidence_stack(lambda,1,1,0,pol,dbr);
            near(std::abs(parity.reflection-br.reflection),0,1e-12,"normal amplitude parity");
            near(std::abs(parity.transmission-br.transmission),0,1e-12,"normal transmission parity");
            const auto tir = oblique_incidence_stack(lambda,1.5,1,1.0,pol,{});
            near(tir.reflectance,1,1e-13,"TIR");
            near(tir.transmittance,0,1e-13,"evanescent exit flux");
            const auto tunnel = oblique_incidence_stack(lambda,1.5,1.5,1.0,pol,gap);
            const auto wide = oblique_incidence_stack(lambda,1.5,1.5,1.0,pol,wide_gap);
            require(tunnel.transmittance>0.1 && tunnel.transmittance<1,"frustrated TIR");
            require(wide.transmittance<tunnel.transmittance*1e-6,"evanescent attenuation");
            near(tunnel.reflectance+tunnel.transmittance,1,1e-12,"tunnel balance");
            const auto slab = oblique_incidence_stack(lambda,1,1,theta,pol,film);
            // Independent symmetric Fabry-Perot intensity formula.
            const double q0=ci, q1=std::sqrt(4-std::pow(std::sin(theta),2));
            const double y0=pol==Polarization::te?q0:1/q0;
            const double y1=pol==Polarization::te?q1:4/q1;
            const double rr=std::pow((y0-y1)/(y0+y1),2);
            const double F=4*rr/std::pow(1-rr,2);
            const double expected=1/(1+F*std::pow(std::sin(k*q1*film[0].thickness_m),2));
            near(slab.transmittance,expected,1e-12,"Fabry-Perot T");
            near(slab.reflectance+slab.transmittance,1,1e-12,"slab balance");
            const auto oblique_loss=oblique_incidence_stack(lambda,1,1,theta,pol,lossy);
            near(oblique_loss.transmittance,0,1e-12,"thick oblique absorption T");
            require(oblique_loss.absorptance>0,"thick oblique absorption");
            const auto skipped=oblique_incidence_stack(lambda,1.5,1.5,1.0,pol,zero);
            near(skipped.transmittance,1,1e-13,"zero thickness layer");
        }
        bool bad_angle=false;
        try { (void)oblique_incidence_stack(lambda,1,1,-1,Polarization::te,{}); }
        catch (const std::invalid_argument&) { bad_angle=true; }
        require(bad_angle,"reject negative angle");
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
