#include "nova/physics/optical_layer_power.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace nova::optics {
namespace {
struct Local {
    Complex propagation;
    Complex right_reflection;
    Complex left_reflection;
};
Complex fresnel(Complex a, Complex b) {
    return (a-b)/(a+b);
}
double flux(Complex forward, Complex backward, Complex index, double incident) {
    const auto e = forward + backward;
    const auto h = index * (forward - backward);
    const double value = std::real(e * std::conj(h)) / incident;
    if (!std::isfinite(value)) throw std::overflow_error("nonfinite layer power");
    return value;
}
}

LayerPowerAudit audit_layer_power(double wavelength, double incident, double exit,
                                  std::span<const Layer> layers) {
    // Reuse the public solver's input validation and external power solution.
    LayerPowerAudit result{normal_incidence_stack(wavelength,incident,exit,layers),{},0.0};
    const double k = 2.0 * std::numbers::pi_v<double> / wavelength;
    std::vector<Local> local(layers.size());
    Complex next_index = exit;
    Complex next_r = 0.0;
    for (std::size_t j=layers.size(); j>0; --j) {
        const auto& layer = layers[j-1];
        const auto ri = fresnel(layer.index,next_index);
        const auto rr = (ri+next_r)/(1.0+ri*next_r);
        const auto p = std::exp(Complex{0.0,1.0} * layer.index * (k*layer.thickness_m));
        const auto rl = p*p*rr;
        local[j-1] = {p,rr,rl};
        next_index = layer.index;
        next_r = rl;
    }
    result.layers.reserve(layers.size());
    Complex previous_index = incident;
    Complex forward = 1.0;
    double total_absorption = 0.0;
    for (std::size_t j=0; j<layers.size(); ++j) {
        const auto& layer = layers[j];
        const auto& data = local[j];
        const auto ri = fresnel(previous_index,layer.index);
        const auto f_left = forward*(1.0+ri)/(1.0+ri*data.left_reflection);
        const auto f_right = f_left*data.propagation;
        const double entering = flux(f_left,f_left*data.left_reflection,layer.index,incident);
        const double leaving = flux(f_right,f_right*data.right_reflection,layer.index,incident);
        const double absorbed = entering-leaving;
        result.layers.push_back({entering,leaving,absorbed});
        total_absorption += absorbed;
        forward = f_right;
        previous_index = layer.index;
    }
    result.balance_defect = total_absorption-result.stack.absorptance;
    if (!std::isfinite(result.balance_defect))
        throw std::overflow_error("nonfinite optical power balance");
    return result;
}
} // namespace nova::optics
