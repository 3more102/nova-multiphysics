#include "test_framework.hpp"
#include "nova/physics/optical_dispersion.hpp"
#include "nova/physics/maxwell_3d.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

using namespace nova::physics;

NOVA_TEST(optics_drude_negative_permittivity_and_passive_loss) {
    const std::array<OpticalOscillator,1> oscillators{{{3.0,0.0,0.5}}};
    const auto e=optical_permittivity(1.0,2.0,oscillators);
    NOVA_REQUIRE_NEAR(e.real(),-5.2,1e-12);
    NOVA_REQUIRE_NEAR(e.imag(),3.6,1e-12);
    const auto high=optical_permittivity(1e6,2.0,oscillators);
    NOVA_REQUIRE_NEAR(high.real(),2.0,1e-10);
    NOVA_REQUIRE(high.imag()>0.0);
}
NOVA_TEST(optics_lorentz_resonance_and_superposition) {
    const std::array<OpticalOscillator,1> one{{{3.0,2.0,0.5}}};
    const std::array<OpticalOscillator,2> two{{{3.0,2.0,0.5},{1.0,0.0,0.0}}};
    const auto resonant=optical_permittivity(2.0,1.5,one);
    NOVA_REQUIRE_NEAR(resonant.real(),1.5,1e-12);
    NOVA_REQUIRE_NEAR(resonant.imag(),9.0,1e-12);
    const auto sum=optical_permittivity(2.0,1.5,two);
    NOVA_REQUIRE_NEAR(sum.real(),1.25,1e-12);
    NOVA_REQUIRE_NEAR(sum.imag(),9.0,1e-12);
    const std::array<OpticalOscillator,1> scaled{{{3e150,2e150,5e149}}};
    const auto large=optical_permittivity(2e150,1.5,scaled);
    NOVA_REQUIRE_NEAR(large.real(),resonant.real(),1e-12);
    NOVA_REQUIRE_NEAR(large.imag(),resonant.imag(),1e-12);
}
NOVA_TEST(optics_invalid_inputs_and_exact_poles) {
    const std::array<OpticalOscillator,1> pole{{{1.0,2.0,0.0}}};
    const std::array<OpticalOscillator,1> gain{{{1.0,2.0,-0.1}}};
    bool rejected=false;
    try {(void)optical_permittivity(2.0,1.0,pole);}
    catch(const std::domain_error&){rejected=true;}
    NOVA_REQUIRE(rejected);
    rejected=false;
    try {(void)optical_permittivity(2.0,1.0,gain);}
    catch(const std::invalid_argument&){rejected=true;}
    NOVA_REQUIRE(rejected);
    rejected=false;
    try {(void)optical_permittivity(0.0,1.0,pole);}
    catch(const std::invalid_argument&){rejected=true;}
    NOVA_REQUIRE(rejected);
    rejected=false;
    try {(void)optical_angular_frequency(std::numeric_limits<double>::infinity());}
    catch(const std::invalid_argument&){rejected=true;}
    NOVA_REQUIRE(rejected);
}
NOVA_TEST(optics_wavelength_conversion_and_maxwell_material_bridge) {
    constexpr double lambda=1550e-9,c=299792458.0;
    const auto omega=optical_angular_frequency(lambda);
    NOVA_REQUIRE_NEAR(omega*lambda/(2*std::numbers::pi*c),1.0,1e-15);
    const std::array<OpticalOscillator,1> metal{{{3*omega,0.0,0.1*omega}}};
    const MaxwellCellMaterial3D material{optical_permittivity(omega,1.0,metal),1.0};
    NOVA_REQUIRE(material.relative_permittivity.real()<0.0);
    NOVA_REQUIRE(material.relative_permittivity.imag()>0.0);
    const auto mesh=nova::fem::make_box_mesh_3d(0,lambda,0,lambda,0,lambda,2,2,2).mesh;
    const std::vector<MaxwellCellMaterial3D> materials(mesh.elements.size(),material);
    const auto system=assemble_maxwell_3d(mesh,materials,omega/c,true);
    NOVA_REQUIRE(!system.free_edges.empty());
    NOVA_REQUIRE(system.real_block.rows()==2*system.free_edges.size());
}
