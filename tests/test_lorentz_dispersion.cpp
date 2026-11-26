#include "test_framework.hpp"
#include "nova/physics/lorentz_dispersion.hpp"
#include <array>
#include <limits>

using namespace nova::optics;

NOVA_TEST(lorentz_static_resonance_and_high_frequency_limits) {
    const std::array poles{LorentzPole{2.0,0.5,3.0}};
    const auto dc = lorentz_permittivity(0.0,2.0,poles);
    NOVA_REQUIRE_NEAR(dc.relative_permittivity.real(),5.0,1e-14);
    NOVA_REQUIRE_NEAR(dc.relative_permittivity.imag(),0.0,1e-14);
    NOVA_REQUIRE_NEAR(dc.derivative_per_angular_frequency.imag(),0.375,1e-14);
    const auto peak = lorentz_permittivity(2.0,2.0,poles);
    NOVA_REQUIRE_NEAR(peak.relative_permittivity.real(),2.0,1e-14);
    NOVA_REQUIRE_NEAR(peak.relative_permittivity.imag(),12.0,1e-13);
    NOVA_REQUIRE_NEAR(peak.derivative_per_angular_frequency.real(),-48.0,1e-12);
    NOVA_REQUIRE_NEAR(peak.derivative_per_angular_frequency.imag(),-6.0,1e-12);
    const auto high = lorentz_permittivity(1e8,2.0,poles);
    NOVA_REQUIRE_NEAR(high.relative_permittivity.real(),2.0,1e-12);
    NOVA_REQUIRE(high.relative_permittivity.imag()>=0.0);
}
NOVA_TEST(lorentz_multipole_passivity_and_derivative) {
    const std::array poles{LorentzPole{2e15,1e14,1.2},LorentzPole{4e15,3e14,0.8}};
    constexpr double h=1e9;
    for (double w : {1e15,2e15,3e15,4e15,6e15}) {
        const auto value=lorentz_permittivity(w,1.5,poles);
        const auto left=lorentz_permittivity(w-h,1.5,poles);
        const auto right=lorentz_permittivity(w+h,1.5,poles);
        const auto numerical=(right.relative_permittivity-left.relative_permittivity)/(2*h);
        NOVA_REQUIRE(value.relative_permittivity.imag()>=0.0);
        NOVA_REQUIRE(std::abs(numerical-value.derivative_per_angular_frequency) <=
            1e-7*std::abs(value.derivative_per_angular_frequency)+1e-24);
    }
}
NOVA_TEST(lorentz_rejects_nonpassive_and_singular_parameters) {
    for (const auto pole : {LorentzPole{2,-1,1},LorentzPole{2,1,-1},
            LorentzPole{0,1,1},LorentzPole{2,1,std::numeric_limits<double>::infinity()}}) {
        bool rejected=false;
        try { (void)lorentz_permittivity(1,1,std::array{pole}); }
        catch (const std::invalid_argument&) { rejected=true; }
        NOVA_REQUIRE(rejected);
    }
    bool rejected=false;
    try { (void)lorentz_permittivity(2,1,std::array{LorentzPole{2,0,1}}); }
    catch (const std::domain_error&) { rejected=true; }
    NOVA_REQUIRE(rejected);
    NOVA_REQUIRE_NEAR(lorentz_permittivity(2,1,std::array{LorentzPole{2,0,0}})
        .relative_permittivity.real(),1.0,0.0);
}
