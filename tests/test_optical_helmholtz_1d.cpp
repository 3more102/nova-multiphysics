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


void exact_zero_permittivity_layer() {
    // At epsilon_r=0 the Helmholtz equation becomes E''=0.
    // The exact affine field is in the P1 trial space on every mesh.
    constexpr double ni=1.3, nt=2.1, length=0.41;
    const Complex I{0.0,1.0};
    const double L=2.0*std::numbers::pi_v<double>*length;
    const Complex t=2.0/(1.0+nt/ni-I*nt*L);
    const Complex r=t*(1.0-I*nt*L)-1.0;
    for (std::size_t count : {std::size_t{1},std::size_t{7},std::size_t{31}}) {
        std::vector<double> nodes(count+1);
        for (std::size_t j=0;j<=count;++j) {
            const double u=static_cast<double>(j)/static_cast<double>(count);
            nodes[j]=wavelength*length*u*u;
        }
        const std::vector<Complex> eps(count,Complex{});
        const auto fem=solve_helmholtz_1d(wavelength,ni,nt,nodes,eps);
        require(std::abs(fem.reflection-r)<1e-10,"ENZ exact complex reflection");
        require(std::abs(fem.transmission-t)<1e-10,"ENZ exact complex transmission");
        for (std::size_t j=0;j<=count;++j) {
            const double s=2.0*std::numbers::pi_v<double>*nodes[j]/wavelength;
            require(std::abs(fem.electric_field[j]-t*(1.0+I*nt*(s-L)))<1e-10,
                    "ENZ exact affine electric field");
        }
        require(fem.absorptance==0.0,"ENZ lossless material");
        for (double a : fem.cell_absorptance) require(a==0.0,"ENZ zero cell loss");
        require(std::abs(fem.reflectance-std::norm(r))<1e-10,"ENZ reflected power");
        require(std::abs(fem.transmittance-(nt/ni)*std::norm(t))<1e-10,
                "ENZ transmitted power");
        require(std::abs(fem.power_balance_residual)<1e-10,"ENZ power balance");
    }
}

void graded_multilayer_reference() {
    // Independent exact propagation of [E,dE/ds], s=k0*x, rather than
    // the scattering cascade used by the existing single-slab tests.
    constexpr double ni=1.3, nt=1.8;
    constexpr double pi=std::numbers::pi_v<double>;
    const Complex I{0,1};
    const std::array<Complex,3> index{{1.55,{0.25,2.0},{2.2,0.08}}};
    const std::array<double,3> width{{0.18,0.07,0.31}};
    std::array<Complex,4> exact_e{}, exact_d{};
    exact_e.back()=1.0;
    exact_d.back()=I*nt;
    for (std::size_t j=index.size(); j-- > 0;) {
        const auto n=index[j];
        const auto phase=n*(2*pi*width[j]);
        const auto c=std::cos(phase), s=std::sin(phase);
        exact_e[j]=c*exact_e[j+1]-s*exact_d[j+1]/n;
        exact_d[j]=n*s*exact_e[j+1]+c*exact_d[j+1];
    }
    const Complex exact_t=2.0/(exact_e.front()+exact_d.front()/(I*ni));
    for (std::size_t j=0;j<exact_e.size();++j) {
        exact_e[j]*=exact_t;
        exact_d[j]*=exact_t;
    }
    const Complex exact_r=exact_e.front()-1.0;
    const auto flux=[&](std::size_t j) {
        return std::imag(std::conj(exact_e[j])*exact_d[j])/ni;
    };
    const auto run=[&](std::size_t m, double origin, bool reverse) {
        std::vector<double> nodes{origin};
        std::vector<Complex> eps;
        for(std::size_t j=0;j<index.size();++j) {
            const std::size_t layer=reverse ? index.size()-1-j : j;
            const double start=nodes.back();
            for(std::size_t q=1;q<=m;++q) {
                const double u=static_cast<double>(q)/static_cast<double>(m);
                nodes.push_back(start+wavelength*width[layer]*std::pow(u,1.25));
                eps.push_back(index[layer]*index[layer]);
            }
        }
        return solve_helmholtz_1d(wavelength,reverse?nt:ni,reverse?ni:nt,nodes,eps);
    };
    const auto coarse=run(80,0,false);
    const auto fine=run(160,0,false);
    require(std::abs(fine.transmission-exact_t)<0.30*std::abs(coarse.transmission-exact_t),
            "graded multilayer second-order transmission convergence");
    require(std::abs(fine.reflection-exact_r)<3e-5,"graded multilayer complex reflection");
    require(std::abs(fine.transmission-exact_t)<1e-4,"graded multilayer complex transmission");
    for (std::size_t j=0;j<exact_e.size();++j)
        require(std::abs(fine.electric_field[j*160]-exact_e[j])<1e-4,
                "graded multilayer interface field");
    for (std::size_t j=0;j<index.size();++j) {
        double loss=0;
        for(std::size_t q=0;q<160;++q) {
            const double cell=fine.cell_absorptance[j*160+q];
            require(cell>=0,"graded multilayer passive cell");
            loss+=cell;
        }
        require(std::abs(loss-(flux(j)-flux(j+1)))<1e-5,
                "layer absorption matches exact flux drop");
    }
    require(std::abs(fine.power_balance_residual)<1e-10,"graded multilayer energy balance");
    const auto translated=run(160,10*wavelength,false);
    require(std::abs(translated.transmission-fine.transmission)<1e-10,
            "mesh translation preserves transmission reference planes");
    const auto reversed=run(160,0,true);
    // Reciprocity is exact for mirrored meshes; grading here is rebuilt
    // from each layer's new left endpoint, so compare at discretization accuracy.
    require(std::abs(reversed.transmittance-fine.transmittance)<2e-4,
            "asymmetric reciprocal multilayer power transmission");
}

}
int main() {
    try {
        graded_multilayer_reference();
        exact_zero_permittivity_layer();
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
