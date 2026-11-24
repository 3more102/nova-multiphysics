#include "test_framework.hpp"
#include "nova/physics/slab_waveguide.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

using namespace nova::physics;
namespace {
SlabWaveguide at_v(double v) {
    return {1.0,v/(2.0*std::numbers::pi_v<double>*std::sqrt(1.25)),1.5,1.0};
}
}

NOVA_TEST(slab_te_tm_independent_reference) {
    const auto g=at_v(2.5);
    const auto te=solve_symmetric_slab(g,SlabPolarization::te);
    const auto tm=solve_symmetric_slab(g,SlabPolarization::tm);
    NOVA_REQUIRE(te.size()==2 && tm.size()==2);
    NOVA_REQUIRE_NEAR(te[0].effective_index,1.4153986308765596,1e-11);
    NOVA_REQUIRE_NEAR(te[1].effective_index,1.160423010582448,1e-11);
    NOVA_REQUIRE_NEAR(tm[0].effective_index,1.3817754032255871,1e-11);
    NOVA_REQUIRE_NEAR(tm[1].effective_index,1.085752934406741,1e-11);
    for (const auto& modes : {te,tm}) {
        for (const auto& mode : modes) {
            const double u=mode.normalized_transverse, w=mode.normalized_decay;
            NOVA_REQUIRE_NEAR(u*u+w*w,6.25,1e-12);
            NOVA_REQUIRE_NEAR(mode.dispersion_residual,0.0,1e-12);
            NOVA_REQUIRE(mode.core_power_fraction>0 && mode.core_power_fraction<1);
            NOVA_REQUIRE_NEAR(mode.propagation_constant,
                2*std::numbers::pi_v<double>*mode.effective_index,1e-12);
        }
        NOVA_REQUIRE(modes[0].core_power_fraction>modes[1].core_power_fraction);
    }
}

NOVA_TEST(slab_mode_cutoff_and_scale) {
    for (auto pol : {SlabPolarization::te,SlabPolarization::tm}) {
        NOVA_REQUIRE(solve_symmetric_slab(at_v(1.56),pol).size()==1);
        NOVA_REQUIRE(solve_symmetric_slab(at_v(1.58),pol).size()==2);
        const auto weak=solve_symmetric_slab(at_v(0.1),pol);
        const auto strong=solve_symmetric_slab(at_v(5.0),pol);
        NOVA_REQUIRE(weak[0].core_power_fraction<strong[0].core_power_fraction);
        auto g=at_v(2.5);
        const auto a=solve_symmetric_slab(g,pol);
        g.vacuum_wavelength*=1e-6; g.half_thickness*=1e-6;
        const auto b=solve_symmetric_slab(g,pol);
        NOVA_REQUIRE_NEAR(a[0].effective_index,b[0].effective_index,1e-13);
        NOVA_REQUIRE_NEAR(a[0].propagation_constant,b[0].propagation_constant*1e-6,1e-11);
    }
}

NOVA_TEST(slab_rejects_unguiding_and_nonpositive_parameters) {
    auto g=at_v(2.5);
    g.core_index=g.cladding_index;
    bool rejected=false;
    try { (void)solve_symmetric_slab(g,SlabPolarization::te); }
    catch (const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected);
    g=at_v(2.5); g.vacuum_wavelength=0;
    rejected=false;
    try { (void)solve_symmetric_slab(g,SlabPolarization::tm); }
    catch (const std::invalid_argument&) { rejected=true; }
    NOVA_REQUIRE(rejected);
}
