#pragma once
#include "nova/physics/optics.hpp"
#include <vector>

namespace nova::optics {
struct LayerPower {
    double entering_flux;
    double leaving_flux;
    double absorptance;
};
struct LayerPowerAudit {
    StackResult stack;
    std::vector<LayerPower> layers;
    double balance_defect;
};
// Unit incident electric field, normal incidence, passive isotropic nonmagnetic layers.
// Fluxes normalized to incident power; layer order is incident to exit.
// Signed roundoff-level defects are retained, not clipped.
[[nodiscard]] LayerPowerAudit audit_layer_power(double wavelength_m,
    double incident_index, double exit_index, std::span<const Layer> layers);
} // namespace nova::optics
