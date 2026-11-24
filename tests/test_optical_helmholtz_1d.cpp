#include "nova/physics/optical_helmholtz_1d.hpp"
#include "nova/physics/optics.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>
#include <utility>

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

void repeated_pivot_dense_reference() {
    // Coarse, deliberately under-resolved meshes test algebra, not wave accuracy.
    // Assemble a dense matrix by two-point Gauss quadrature and solve it with
    // ordinary full-row partial pivoting, independently of the banded solver.
    constexpr std::size_t n=13;
    constexpr double ni=0.1, nt=0.2;
    constexpr double pi=std::numbers::pi_v<double>;
    const Complex I{0,1};
    const double gauss=1.0/std::sqrt(3.0);
    for (bool lossy : {false,true}) {
        std::vector<double> x(n);
        std::vector<Complex> eps(n-1);
        std::vector<std::vector<Complex>> a(n,std::vector<Complex>(n));
        std::vector<Complex> b(n);
        for (std::size_t j=0;j+1<n;++j) {
            const double step=std::sqrt(3.0)*(lossy ? 1.0+0.03*static_cast<double>(j%3) : 1.0);
            x[j+1]=x[j]+step*wavelength/(2*pi);
            eps[j]=lossy ? Complex{1.0,0.02*static_cast<double>(1+j%4)} : Complex{1.0};
            const double h=(x[j+1]-x[j])*2*pi/wavelength;
            const std::array<double,2> grad{-1.0/h,1.0/h};
            for (double xi : {-gauss,gauss}) {
                const std::array<double,2> phi{(1-xi)/2,(1+xi)/2};
                for (std::size_t p=0;p<2;++p)
                    for (std::size_t q=0;q<2;++q)
                        a[j+p][j+q]+=(h/2)*(grad[p]*grad[q]-eps[j]*phi[p]*phi[q]);
            }
        }
        a.front().front()-=I*ni;
        a.back().back()-=I*nt;
        b.front()=-2.0*I*ni;
        const auto original=a;
        const auto load=b;
        std::size_t swaps=0;
        for (std::size_t k=0;k<n;++k) {
            std::size_t pivot=k;
            for (std::size_t p=k+1;p<n;++p)
                if (std::abs(a[p][k])>std::abs(a[pivot][k])) pivot=p;
            require(std::abs(a[pivot][k])>1e-12,"nonsingular dense reference");
            if (pivot!=k) {
                std::swap(a[k],a[pivot]);
                std::swap(b[k],b[pivot]);
                ++swaps;
            }
            for (std::size_t p=k+1;p<n;++p) {
                const Complex factor=a[p][k]/a[k][k];
                for (std::size_t q=k;q<n;++q) a[p][q]-=factor*a[k][q];
                b[p]-=factor*b[k];
            }
        }
        require(swaps>=2,"reference exercises repeated pivoting");
        for (std::size_t k=n;k-- >0;) {
            for (std::size_t q=k+1;q<n;++q) b[k]-=a[k][q]*b[q];
            b[k]/=a[k][k];
        }
        const auto fem=solve_helmholtz_1d(wavelength,ni,nt,x,eps);
        for (std::size_t p=0;p<n;++p) {
            require(std::abs(fem.electric_field[p]-b[p])<1e-11,
                    "banded field matches independent dense solve");
            Complex residual=-load[p];
            for (std::size_t q=0;q<n;++q) residual+=original[p][q]*fem.electric_field[q];
            require(std::abs(residual)<1e-11,"independent assembled equation residual");
        }
        require(std::abs(fem.reflection-(b.front()-1.0))<1e-11,"dense reference reflection");
        require(std::abs(fem.transmission-b.back())<1e-11,"dense reference transmission");
        require(std::abs(fem.power_balance_residual)<1e-11,"repeated pivot energy balance");
    }
}

}
int main() {
    try {
        graded_multilayer_reference();
        repeated_pivot_dense_reference();
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
