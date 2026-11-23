#include "nova/physics/transient_electrothermal_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace nova::physics {
namespace {

void require_same_mesh(const fem::TetrahedralMesh3D& a,
                       const fem::TetrahedralMesh3D& b) {
    if (a.nodes.size() != b.nodes.size() ||
        a.elements.size() != b.elements.size()) {
        throw std::invalid_argument("Transient electrothermal 3D meshes differ in size.");
    }
    for (std::size_t i = 0; i < a.nodes.size(); ++i) {
        if (a.nodes[i].x != b.nodes[i].x ||
            a.nodes[i].y != b.nodes[i].y ||
            a.nodes[i].z != b.nodes[i].z) {
            throw std::invalid_argument("Transient electrothermal 3D mesh nodes differ.");
        }
    }
    for (std::size_t i = 0; i < a.elements.size(); ++i) {
        if (a.elements[i] != b.elements[i]) {
            throw std::invalid_argument("Transient electrothermal 3D element order differs.");
        }
    }
}

void require_cell_material(const fem::ScalarDiffusion3DProblem& p) {
    const bool scalar = !p.cell_diffusivity.empty();
    const bool tensor = !p.cell_diffusivity_tensors.empty();
    if (scalar == tensor || p.diffusivity || p.diffusivity_tensor ||
        (scalar && p.cell_diffusivity.size() != p.mesh.elements.size()) ||
        (tensor && p.cell_diffusivity_tensors.size() != p.mesh.elements.size())) {
        throw std::invalid_argument(
            "Temperature-dependent transient material requires one complete "
            "per-cell scalar or SPD-tensor reference conductivity.");
    }
}

void update_material(
    const fem::ScalarDiffusion3DProblem& reference,
    fem::ScalarDiffusion3DProblem& target,
    const std::vector<double>& temperatures,
    const std::function<double(double)>& scale_at_temperature
) {
    if (!scale_at_temperature) return;
    for (std::size_t i = 0; i < reference.mesh.elements.size(); ++i) {
        const auto& tetra = reference.mesh.elements[i];
        long double sum = 0.0L;
        for (const std::size_t node : tetra) {
            if (node >= temperatures.size())
                throw std::invalid_argument("Transient electrothermal node is out of range.");
            sum += static_cast<long double>(temperatures[node]);
        }
        const double mean = static_cast<double>(sum / 4.0L);
        if (!std::isfinite(mean))
            throw std::overflow_error("Transient electrothermal cell mean overflows.");
        const double factor = scale_at_temperature(mean);
        if (!std::isfinite(factor) || !(factor > 0.0))
            throw std::invalid_argument("Transient electrothermal material scale must be positive.");
        const auto multiply = [factor](double coefficient) {
            const double value = coefficient * factor;
            if (!std::isfinite(value))
                throw std::overflow_error("Transient electrothermal material overflows.");
            return value;
        };
        if (!reference.cell_diffusivity.empty()) {
            const double value = multiply(reference.cell_diffusivity[i]);
            if (!(value > 0.0))
                throw std::invalid_argument("Transient electrothermal scalar material underflows.");
            target.cell_diffusivity[i] = value;
        } else {
            const auto& a = reference.cell_diffusivity_tensors[i];
            target.cell_diffusivity_tensors[i] = {
                multiply(a.xx), multiply(a.yy), multiply(a.zz),
                multiply(a.xy), multiply(a.xz), multiply(a.yz)
            };
            // The underlying FEM solver performs the final SPD validation.
        }
    }
}

void add_joule_heat(std::vector<double>& source, const std::vector<double>& joule) {
    if (!source.empty() && source.size() != joule.size()) {
        throw std::invalid_argument("Transient electrothermal cell-source length mismatch.");
    }
    if (source.empty()) source.assign(joule.size(), 0.0);
    for (std::size_t i = 0; i < joule.size(); ++i) {
        if (!std::isfinite(source[i]) || !std::isfinite(joule[i])) {
            throw std::invalid_argument("Transient electrothermal source is not finite.");
        }
        const long double combined = static_cast<long double>(source[i]) + joule[i];
        if (!std::isfinite(combined) ||
            std::abs(combined) >
                static_cast<long double>(std::numeric_limits<double>::max())) {
            throw std::overflow_error("Transient electrothermal heat source overflows.");
        }
        source[i] = static_cast<double>(combined);
    }
}

void validate(const TransientElectrothermal3DProblem& p,
              const TransientElectrothermal3DOptions& o) {
    require_same_mesh(p.electrical.mesh, p.thermal.diffusion.mesh);
    if (p.electrical_conductivity_scale) require_cell_material(p.electrical);
    if (p.thermal_conductivity_scale) {
        require_cell_material(p.thermal.diffusion);
        if (p.thermal.time_dependent_diffusivity ||
            p.thermal.time_dependent_diffusivity_tensor) {
            throw std::invalid_argument(
                "Thermal k(T) cannot be combined with time-varying thermal conductivity.");
        }
    }
    if (o.thermal_time.scheme != fem::TransientHeat3DScheme::BackwardEuler ||
        o.thermal_time.on_step ||
        o.thermal_time.steps == 0 ||
        !std::isfinite(o.thermal_time.start_time) ||
        !std::isfinite(o.thermal_time.time_step) ||
        !(o.thermal_time.time_step > 0.0)) {
        throw std::invalid_argument(
            "Transient electrothermal coupling requires fixed-step BE and no on_step callback.");
    }
    if (o.max_coupling_iterations == 0 ||
        !std::isfinite(o.relaxation) || !(o.relaxation > 0.0) ||
        o.relaxation > 1.0 ||
        !std::isfinite(o.absolute_tolerance_kelvin) ||
        o.absolute_tolerance_kelvin < 0.0 ||
        !std::isfinite(o.relative_tolerance) ||
        o.relative_tolerance < 0.0 ||
        (o.absolute_tolerance_kelvin == 0.0 && o.relative_tolerance == 0.0)) {
        throw std::invalid_argument("Invalid transient electrothermal Picard parameters.");
    }
    if (p.thermal.initial_temperature.size() != p.electrical.mesh.nodes.size()) {
        throw std::invalid_argument("Transient electrothermal initial state length mismatch.");
    }
    for (const auto t : p.thermal.initial_temperature) {
        if (!std::isfinite(t))
            throw std::invalid_argument("Transient electrothermal initial state is nonfinite.");
    }
    const double final = o.thermal_time.start_time +
        o.thermal_time.time_step * static_cast<double>(o.thermal_time.steps);
    if (!std::isfinite(final) || !(final > o.thermal_time.start_time)) {
        throw std::invalid_argument("Transient electrothermal time interval is invalid.");
    }
}

} // namespace

TransientElectrothermal3DResult solve_transient_electrothermal_3d(
    const TransientElectrothermal3DProblem& p,
    const TransientElectrothermal3DOptions& o
) {
    validate(p, o);
    TransientElectrothermal3DResult result;
    result.temperature = p.thermal.initial_temperature;
    result.final_time = o.thermal_time.start_time;
    result.accepted_steps.reserve(o.thermal_time.steps);

    for (std::size_t step = 0; step < o.thermal_time.steps; ++step) {
        const double begin = o.thermal_time.start_time +
            o.thermal_time.time_step * static_cast<double>(step);
        const double end = o.thermal_time.start_time +
            o.thermal_time.time_step * static_cast<double>(step + 1U);
        if (!std::isfinite(begin) || !std::isfinite(end) || !(end > begin)) {
            throw std::invalid_argument("Transient electrothermal step does not advance time.");
        }

        std::vector<double> iterate = result.temperature;
        double last_defect = 0.0;
        double last_normalized = 0.0;
        bool accepted = false;
        for (std::size_t iteration = 1; iteration <= o.max_coupling_iterations; ++iteration) {
            auto electrical = p.electrical;
            update_material(p.electrical, electrical, iterate,
                            p.electrical_conductivity_scale);
            auto joule = solve_dc_joule_heating_3d(electrical,
                                                   o.electrical_solver_options);

            auto thermal = p.thermal;
            thermal.initial_temperature = result.temperature;  // T^n, not the trial
            update_material(p.thermal.diffusion, thermal.diffusion, iterate,
                            p.thermal_conductivity_scale);
            add_joule_heat(thermal.diffusion.cell_source, joule.joule_density);

            auto heat_options = o.thermal_time;
            heat_options.start_time = begin;
            heat_options.steps = 1U;
            heat_options.retain_step_solutions = false;
            auto heat = fem::solve_transient_heat_3d(thermal, heat_options);
            if (heat.linear_solvers.size() != 1U ||
                !heat.linear_solvers.front().converged ||
                heat.diagnostics.size() != 1U ||
                heat.temperature.size() != iterate.size()) {
                throw std::runtime_error("Transient electrothermal thermal solve failed.");
            }
            last_defect = 0.0;
            last_normalized = 0.0;
            for (std::size_t node = 0; node < iterate.size(); ++node) {
                const double next = heat.temperature[node];
                const double old = iterate[node];
                if (!std::isfinite(next))
                    throw std::runtime_error("Transient electrothermal temperature is nonfinite.");
                const double defect = std::abs(next - old);
                const double norm = o.absolute_tolerance_kelvin +
                    o.relative_tolerance * std::max(std::abs(next), std::abs(old));
                if (!std::isfinite(norm) || !(norm > 0.0) ||
                    !std::isfinite(defect / norm)) {
                    throw std::overflow_error("Transient electrothermal residual overflows.");
                }
                last_defect = std::max(last_defect, defect);
                last_normalized = std::max(last_normalized, defect / norm);
            }
            if (last_normalized <= 1.0) {
                result.temperature = std::move(heat.temperature);
                result.final_time = end;
                result.last_electrical = std::move(joule);
                result.accepted_steps.push_back({
                    end, iteration, last_defect, last_normalized,
                    result.last_electrical.total_joule_power_watts, heat.diagnostics.front()
                });
                accepted = true;
                break;
            }
            if (iteration == o.max_coupling_iterations) break;
            for (std::size_t node = 0; node < iterate.size(); ++node) {
                iterate[node] += o.relaxation * (heat.temperature[node] - iterate[node]);
                if (!std::isfinite(iterate[node])) {
                    throw std::overflow_error("Transient electrothermal relaxation overflows.");
                }
            }
        }
        if (!accepted) {
            result.converged = false;
            result.failed_step = step + 1U;
            result.failed_max_normalized_residual = last_normalized;
            return result;
        }
    }
    return result;
}

} // namespace nova::physics
