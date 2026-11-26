#include "nova/physics/optical_layer_power.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {
void near(double a,double b,double tol,const char* message) {
    if (!std::isfinite(a) || std::abs(a-b)>tol) throw std::runtime_error(message);
}
}
int main() {
    using namespace nova::optics;
    constexpr double wavelength=1.55e-6;
    try {
        const auto empty=audit_layer_power(wavelength,1,1.5,{});
        near(empty.balance_defect,0,1e-12,"empty stack balance");
        const std::array<Layer,2> transparent{{{1.5,wavelength/6},{2.0,wavelength/8}}};
        const auto lossless=audit_layer_power(wavelength,1,1,transparent);
        for (const auto& layer:lossless.layers)
            near(layer.absorptance,0,1e-12,"lossless layer absorption");
        near(lossless.balance_defect,0,1e-12,"lossless balance");
        const std::array<Layer,3> mixed{{{{1.5,0.2},wavelength/3},
            {2.0,wavelength/8},{{0.3,4.0},wavelength/10}}};
        const auto audit=audit_layer_power(wavelength,1,1.2,mixed);
        near(audit.layers.front().entering_flux,1-audit.stack.reflectance,1e-12,"entrance flux");
        near(audit.layers.back().leaving_flux,audit.stack.transmittance,1e-12,"exit flux");
        for (std::size_t j=1;j<audit.layers.size();++j)
            near(audit.layers[j-1].leaving_flux,audit.layers[j].entering_flux,1e-12,"interface continuity");
        if (!(audit.layers[0].absorptance>0 && audit.layers[2].absorptance>0))
            throw std::runtime_error("passive absorbing layers");
        near(audit.layers[1].absorptance,0,1e-12,"transparent middle layer");
        near(audit.balance_defect,0,1e-12,"mixed balance");
        const std::array<Layer,1> thick{{{{1.5,0.5},1e-3}}};
        const auto opaque=audit_layer_power(wavelength,1,1,thick);
        near(opaque.layers[0].leaving_flux,0,1e-12,"opaque transmission");
        near(opaque.layers[0].absorptance,1-opaque.stack.reflectance,1e-12,"opaque absorption");
        const std::array<Layer,1> zero{{{{0.3,4.0},0.0}}};
        const auto vanished=audit_layer_power(wavelength,1,1,zero);
        near(vanished.layers[0].absorptance,0,1e-12,"zero thickness");
        near(vanished.stack.reflectance,0,1e-12,"zero thickness reflection");
        // Splitting a homogeneous film must preserve total power and interface flux.
        const std::array<Layer,1> one{{{{1.7,0.1},wavelength/2}}};
        const std::array<Layer,2> split{{{{1.7,0.1},wavelength/4},{{1.7,0.1},wavelength/4}}};
        const auto a=audit_layer_power(wavelength,1,1,one);
        const auto b=audit_layer_power(wavelength,1,1,split);
        near(a.layers[0].absorptance,b.layers[0].absorptance+b.layers[1].absorptance,1e-12,"split invariance");
        near(b.layers[0].leaving_flux,b.layers[1].entering_flux,1e-12,"split interface");
        near(b.balance_defect,0,1e-12,"split balance");
        std::cout<<"NOVA layer power reference cases passed\n";
    } catch (const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
