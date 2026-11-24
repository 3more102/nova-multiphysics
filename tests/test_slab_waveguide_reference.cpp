#include "nova/physics/slab_waveguide.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace {
void require(bool ok, const char* msg) {
    if (!ok) throw std::runtime_error(msg);
}
void near(double a, double b, double tol, const char* msg) {
    require(std::isfinite(a) && std::abs(a-b)<=tol,msg);
}
}
int main() {
    using namespace nova::optics;
    constexpr double pi=std::numbers::pi_v<double>, lambda=1.55e-6;
    constexpr double nc=2.0, ns=1.5;
    const double k=2*pi/lambda, contrast=nc*nc-ns*ns;
    const auto thickness=[&](double v){return 2*v/(k*std::sqrt(contrast));};
    try {
        // Exact TE0 case u=w=pi/4.
        const auto exact=symmetric_slab_modes(lambda,thickness(pi/std::sqrt(8.0)),nc,ns,SlabPolarization::TE);
        require(exact.size()==1,"single mode");
        near(exact[0].effective_index,std::sqrt((nc*nc+ns*ns)/2),1e-13,"exact TE0");
        near(exact[0].core_transverse_wavenumber_per_m/k,exact[0].cladding_decay_per_m/k,1e-13,"equal transverse constants");
        // Independent tangent-form boundary equations for both parity classes.
        for (auto p : {SlabPolarization::TE,SlabPolarization::TM}) {
            const double d=thickness(5.0);
            const auto modes=symmetric_slab_modes(lambda,d,nc,ns,p);
            require(modes.size()==4,"multimode count");
            double previous=nc;
            for (const auto& m:modes) {
                require(m.effective_index<previous && m.effective_index>ns,"ordered bound modes");
                previous=m.effective_index;
                const double u=m.core_transverse_wavenumber_per_m*d/2;
                const double w=m.cladding_decay_per_m*d/2;
                near(u*u+w*w,25,1e-12,"normalized circle");
                near(m.effective_index*m.effective_index+
                     std::pow(m.core_transverse_wavenumber_per_m/k,2),nc*nc,1e-12,"core dispersion");
                const double q=p==SlabPolarization::TE?1.0:nc*nc/(ns*ns);
                const double lhs=m.order%2==0?u*std::tan(u):-u/std::tan(u);
                near(lhs,q*w,1e-10,"interface continuity");
                require(m.normalized_dispersion_residual<1e-12,"root residual");
            }
        }
        const auto te=symmetric_slab_modes(lambda,thickness(2.0),nc,ns,SlabPolarization::TE);
        const auto tm=symmetric_slab_modes(lambda,thickness(2.0),nc,ns,SlabPolarization::TM);
        require(te[0].effective_index>tm[0].effective_index,"TE/TM ordering");
        require(symmetric_slab_modes(lambda,thickness(pi/2-1e-3),nc,ns,SlabPolarization::TE).size()==1,"below cutoff");
        require(symmetric_slab_modes(lambda,thickness(pi/2+1e-3),nc,ns,SlabPolarization::TE).size()==2,"above cutoff");
        const auto thick=symmetric_slab_modes(lambda,thickness(8.0),nc,ns,SlabPolarization::TE);
        require(thick[0].effective_index>te[0].effective_index,"thickness confinement");
        bool rejected=false;
        try { (void)symmetric_slab_modes(lambda,1e-6,ns,nc,SlabPolarization::TE); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"invalid index contrast");
        rejected=false;
        try { (void)symmetric_slab_modes(std::numeric_limits<double>::quiet_NaN(),1e-6,nc,ns,SlabPolarization::TE); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"nonfinite wavelength");
        std::cout<<"NOVA symmetric slab reference checks passed\n";
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
