#include "nova/io/adaptive_electrothermal_checkpoint_3d.hpp"

#include "nova/io/thermal_checkpoint_3d.hpp"

#include <charconv>
#include <cmath>
#include <iomanip>
#include <initializer_list>
#include <istream>
#include <limits>
#include <locale>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>

namespace nova::io {
namespace {

void expect(std::istream& input, const char* label) {
    std::string token;
    if (!(input >> token) || token != label) {
        throw std::invalid_argument(
            std::string("Invalid adaptive electrothermal checkpoint: expected ") + label);
    }
}

template<typename T>
T read_number(std::istream& input, const char* label) {
    expect(input, label);
    std::string token;
    if (!(input >> token)) {
        throw std::invalid_argument(
            std::string("Missing adaptive electrothermal checkpoint value: ") + label);
    }
    T parsed{};
    const auto r = std::from_chars(token.data(), token.data() + token.size(), parsed);
    if (r.ec != std::errc{} || r.ptr != token.data() + token.size()) {
        throw std::invalid_argument(
            std::string("Invalid adaptive electrothermal checkpoint value: ") + label);
    }
    return parsed;
}

void check_finite(double value, const char* what) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string("Nonfinite adaptive checkpoint ") + what);
    }
}

bool valid_model_id(const std::string& id) {
    if (id.empty() || id.size() > 128U) return false;
    for (const unsigned char c : id) {
        if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
              (c >= '0' && c <= '9') || c == '_' || c == '-' ||
              c == '.' || c == ':')) return false;
    }
    return true;
}

void validate_options(const physics::AdaptiveTransientElectrothermal3DOptions& o) {
    const auto positive = [](double value) {
        return std::isfinite(value) && value > 0.0;
    };
    const auto nonnegative = [](double value) {
        return std::isfinite(value) && value >= 0.0;
    };
    const auto& c = o.coupling;
    const auto& thermal = c.thermal_time;
    const auto& e = c.electrical_solver_options;
    const auto& h = thermal.solver_options;
    const auto valid_preconditioner = [](linalg::PreconditionerType p) {
        return p == linalg::PreconditionerType::None ||
               p == linalg::PreconditionerType::Jacobi ||
               p == linalg::PreconditionerType::IncompleteCholesky0;
    };
    if (!std::isfinite(thermal.start_time) || !std::isfinite(o.end_time) ||
        !(o.end_time > thermal.start_time) ||
        !positive(o.initial_time_step) || !positive(o.minimum_time_step) ||
        !positive(o.maximum_time_step) ||
        o.minimum_time_step > o.initial_time_step ||
        o.initial_time_step > o.maximum_time_step ||
        !nonnegative(o.absolute_error_kelvin) || !nonnegative(o.relative_error) ||
        (o.absolute_error_kelvin == 0.0 && o.relative_error == 0.0) ||
        !positive(o.safety_factor) || o.safety_factor > 1.0 ||
        !positive(o.minimum_shrink_factor) || o.minimum_shrink_factor >= 1.0 ||
        !positive(o.maximum_growth_factor) || o.maximum_growth_factor <= 1.0 ||
        o.maximum_attempts == 0U || o.maximum_accepted_steps == 0U ||
        c.max_coupling_iterations == 0U || !positive(c.relaxation) ||
        c.relaxation > 1.0 ||
        !nonnegative(c.absolute_tolerance_kelvin) ||
        !nonnegative(c.relative_tolerance) ||
        (c.absolute_tolerance_kelvin == 0.0 && c.relative_tolerance == 0.0) ||
        h.max_iterations == 0U || e.max_iterations == 0U ||
        !nonnegative(h.relative_tolerance) ||
        !nonnegative(h.absolute_tolerance) ||
        !nonnegative(e.relative_tolerance) ||
        !nonnegative(e.absolute_tolerance) ||
        !valid_preconditioner(h.preconditioner) ||
        !valid_preconditioner(e.preconditioner) ||
        thermal.scheme != fem::TransientHeat3DScheme::BackwardEuler ||
        static_cast<bool>(thermal.on_step)) {
        throw std::invalid_argument("Invalid adaptive electrothermal checkpoint settings.");
    }
}

template<typename T>
void write_setting(std::ostream& output, const char* key, T value) {
    output << key << ' ' << value << '\n';
}

template<typename T>
void read_matching_setting(std::istream& input, const char* key, T expected) {
    const T saved = read_number<T>(input, key);
    if constexpr (std::is_floating_point_v<T>) {
        check_finite(saved, key);
    }
    if (saved != expected) {
        throw std::invalid_argument(
            std::string("Adaptive electrothermal checkpoint setting mismatch: ") + key);
    }
}

void write_options(std::ostream& out,
                   const physics::AdaptiveTransientElectrothermal3DOptions& o) {
    const auto& c = o.coupling;
    const auto& h = c.thermal_time.solver_options;
    const auto& e = c.electrical_solver_options;
    write_setting(out, "end_time", o.end_time);
    write_setting(out, "min_dt", o.minimum_time_step);
    write_setting(out, "max_dt", o.maximum_time_step);
    write_setting(out, "absolute_error_kelvin", o.absolute_error_kelvin);
    write_setting(out, "relative_error", o.relative_error);
    write_setting(out, "safety_factor", o.safety_factor);
    write_setting(out, "minimum_shrink_factor", o.minimum_shrink_factor);
    write_setting(out, "maximum_growth_factor", o.maximum_growth_factor);
    write_setting(out, "maximum_attempts", o.maximum_attempts);
    write_setting(out, "maximum_accepted_steps", o.maximum_accepted_steps);
    write_setting(out, "coupling_max_iterations", c.max_coupling_iterations);
    write_setting(out, "coupling_relaxation", c.relaxation);
    write_setting(out, "coupling_abs_kelvin", c.absolute_tolerance_kelvin);
    write_setting(out, "coupling_rel_tolerance", c.relative_tolerance);
    write_setting(out, "thermal_max_iterations", h.max_iterations);
    write_setting(out, "thermal_relative_tolerance", h.relative_tolerance);
    write_setting(out, "thermal_absolute_tolerance", h.absolute_tolerance);
    write_setting(out, "thermal_preconditioner", static_cast<int>(h.preconditioner));
    write_setting(out, "electrical_max_iterations", e.max_iterations);
    write_setting(out, "electrical_relative_tolerance", e.relative_tolerance);
    write_setting(out, "electrical_absolute_tolerance", e.absolute_tolerance);
    write_setting(out, "electrical_preconditioner", static_cast<int>(e.preconditioner));
}

void read_matching_options(std::istream& in,
                           const physics::AdaptiveTransientElectrothermal3DOptions& o) {
    const auto& c = o.coupling;
    const auto& h = c.thermal_time.solver_options;
    const auto& e = c.electrical_solver_options;
    read_matching_setting(in, "end_time", o.end_time);
    read_matching_setting(in, "min_dt", o.minimum_time_step);
    read_matching_setting(in, "max_dt", o.maximum_time_step);
    read_matching_setting(in, "absolute_error_kelvin", o.absolute_error_kelvin);
    read_matching_setting(in, "relative_error", o.relative_error);
    read_matching_setting(in, "safety_factor", o.safety_factor);
    read_matching_setting(in, "minimum_shrink_factor", o.minimum_shrink_factor);
    read_matching_setting(in, "maximum_growth_factor", o.maximum_growth_factor);
    read_matching_setting(in, "maximum_attempts", o.maximum_attempts);
    read_matching_setting(in, "maximum_accepted_steps", o.maximum_accepted_steps);
    read_matching_setting(in, "coupling_max_iterations", c.max_coupling_iterations);
    read_matching_setting(in, "coupling_relaxation", c.relaxation);
    read_matching_setting(in, "coupling_abs_kelvin", c.absolute_tolerance_kelvin);
    read_matching_setting(in, "coupling_rel_tolerance", c.relative_tolerance);
    read_matching_setting(in, "thermal_max_iterations", h.max_iterations);
    read_matching_setting(in, "thermal_relative_tolerance", h.relative_tolerance);
    read_matching_setting(in, "thermal_absolute_tolerance", h.absolute_tolerance);
    read_matching_setting(in, "thermal_preconditioner", static_cast<int>(h.preconditioner));
    read_matching_setting(in, "electrical_max_iterations", e.max_iterations);
    read_matching_setting(in, "electrical_relative_tolerance", e.relative_tolerance);
    read_matching_setting(in, "electrical_absolute_tolerance", e.absolute_tolerance);
    read_matching_setting(in, "electrical_preconditioner", static_cast<int>(e.preconditioner));
}

void write_progress(
    std::ostream& output,
    const physics::AdaptiveTransientElectrothermal3DProgress& progress
) {
    write_setting(output, "progress_attempted", progress.attempted_steps);
    write_setting(output, "progress_accepted", progress.accepted_steps);
    write_setting(output, "progress_rejected", progress.rejected_steps);
    write_setting(output, "progress_iterations", progress.nonlinear_iterations);
    write_setting(output, "progress_source_energy", progress.source_energy);
    write_setting(output, "progress_neumann_energy", progress.neumann_energy);
    write_setting(output, "progress_robin_energy", progress.robin_energy);
    write_setting(output, "progress_dirichlet_energy", progress.dirichlet_energy);
    write_setting(output, "progress_balance_residual", progress.balance_residual);
    write_setting(output, "progress_heat_content", progress.latest_heat_content);
}

physics::AdaptiveTransientElectrothermal3DProgress
read_progress(std::istream& input) {
    physics::AdaptiveTransientElectrothermal3DProgress progress;
    progress.attempted_steps = read_number<std::size_t>(input, "progress_attempted");
    progress.accepted_steps = read_number<std::size_t>(input, "progress_accepted");
    progress.rejected_steps = read_number<std::size_t>(input, "progress_rejected");
    progress.nonlinear_iterations = read_number<std::size_t>(
        input, "progress_iterations");
    progress.source_energy = read_number<double>(input, "progress_source_energy");
    progress.neumann_energy = read_number<double>(input, "progress_neumann_energy");
    progress.robin_energy = read_number<double>(input, "progress_robin_energy");
    progress.dirichlet_energy = read_number<double>(input, "progress_dirichlet_energy");
    progress.balance_residual = read_number<double>(input, "progress_balance_residual");
    progress.latest_heat_content = read_number<double>(input, "progress_heat_content");
    return progress;
}

void validate_progress(const physics::AdaptiveTransientElectrothermal3DProgress& p) {
    if (p.accepted_steps > p.attempted_steps ||
        p.rejected_steps != p.attempted_steps - p.accepted_steps) {
        throw std::invalid_argument("Invalid adaptive checkpoint cumulative counts.");
    }
    for (const double v : {p.source_energy, p.neumann_energy, p.robin_energy,
                           p.dirichlet_energy, p.balance_residual,
                           p.latest_heat_content}) {
        check_finite(v, "cumulative energy");
    }
}

void validate_state(const physics::AdaptiveTransientElectrothermal3DOptions& o,
                    const AdaptiveElectrothermalCheckpoint3D& state) {
    check_finite(state.time, "time");
    check_finite(state.next_time_step, "next timestep");
    validate_progress(state.progress);
    if (!(state.time >= o.coupling.thermal_time.start_time) ||
        !(state.time <= o.end_time) ||
        state.next_time_step < o.minimum_time_step ||
        state.next_time_step > o.maximum_time_step) {
        throw std::invalid_argument("Adaptive electrothermal checkpoint state out of bounds.");
    }
}

} // namespace

void write_adaptive_electrothermal_checkpoint_3d(
    std::ostream& output,
    const fem::TetrahedralMesh3D& mesh,
    const std::string& model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& options,
    const AdaptiveElectrothermalCheckpoint3D& checkpoint
) {
    if (!valid_model_id(model_id)) {
        throw std::invalid_argument("Invalid adaptive electrothermal checkpoint model identifier.");
    }
    validate_options(options);
    validate_state(options, checkpoint);
    if (!checkpoint.progress_complete) {
        throw std::invalid_argument(
            "Incomplete legacy adaptive checkpoint progress cannot be serialized as v2.");
    }
    if (!output.good()) {
        throw std::runtime_error("Adaptive electrothermal checkpoint output is not writable.");
    }
    std::ostringstream staged;
    staged.imbue(std::locale::classic());
    staged << std::setprecision(std::numeric_limits<double>::max_digits10);
    staged << "NOVA_ADAPTIVE_ELECTROTHERMAL_CHECKPOINT_3D\nversion 2\n";
    staged << "model_id " << model_id << '\n';
    write_options(staged, options);
    write_setting(staged, "next_dt", checkpoint.next_time_step);
    write_progress(staged, checkpoint.progress);
    staged << "thermal_state\n";
    // The pre-existing thermal checkpoint validates the mesh, state size,
    // nonfinite entries and topology, and writes to this temporary stream.
    write_transient_heat_checkpoint_3d(
        staged, mesh, checkpoint.time, checkpoint.temperature);
    if (!staged) {
        throw std::runtime_error("Failed to stage adaptive electrothermal checkpoint.");
    }
    output << staged.str();
    if (!output) {
        throw std::runtime_error("Failed to write adaptive electrothermal checkpoint.");
    }
}

AdaptiveElectrothermalCheckpoint3D read_adaptive_electrothermal_checkpoint_3d(
    std::istream& input,
    const fem::TetrahedralMesh3D& expected_mesh,
    const std::string& expected_model_id,
    const physics::AdaptiveTransientElectrothermal3DOptions& expected_options
) {
    if (!valid_model_id(expected_model_id)) {
        throw std::invalid_argument("Invalid expected adaptive checkpoint model identifier.");
    }
    validate_options(expected_options);
    expect(input, "NOVA_ADAPTIVE_ELECTROTHERMAL_CHECKPOINT_3D");
    const auto version = read_number<unsigned>(input, "version");
    if (version != 1U && version != 2U) {
        throw std::invalid_argument("Unsupported adaptive electrothermal checkpoint version.");
    }
    expect(input, "model_id");
    std::string model_id;
    if (!(input >> model_id) || !valid_model_id(model_id) ||
        model_id != expected_model_id) {
        throw std::invalid_argument("Adaptive electrothermal checkpoint model mismatch.");
    }
    read_matching_options(input, expected_options);
    AdaptiveElectrothermalCheckpoint3D result;
    result.next_time_step = read_number<double>(input, "next_dt");
    if (version == 2U) {
        result.progress = read_progress(input);
    } else {
        result.progress_complete = false;
    }
    expect(input, "thermal_state");
    auto thermal = read_transient_heat_checkpoint_3d(input, expected_mesh);
    result.time = thermal.time;
    result.temperature = std::move(thermal.temperature);
    validate_state(expected_options, result);
    return result;
}

} // namespace nova::io
