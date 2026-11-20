#include "nova/io/thermal_checkpoint_3d.hpp"

#include <charconv>
#include <cmath>
#include <iomanip>
#include <system_error>
#include <istream>
#include <limits>
#include <locale>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace nova::io {
namespace {

void expect(std::istream& input, const char* label) {
    std::string token;
    if (!(input >> token) || token != label) {
        throw std::invalid_argument(std::string("Invalid 3D thermal checkpoint: expected ") + label);
    }
}

template<typename T>
T value(std::istream& input, const char* description) {
    std::string token;
    if (!(input >> token)) {
        throw std::invalid_argument(
            std::string("Invalid 3D thermal checkpoint: cannot read ") + description);
    }
    T parsed{};
    const auto result = std::from_chars(
        token.data(), token.data() + token.size(), parsed);
    if (result.ec != std::errc{} || result.ptr != token.data() + token.size()) {
        throw std::invalid_argument(
            std::string("Invalid 3D thermal checkpoint: cannot read ") + description);
    }
    return parsed;
}

void check_finite(double input, const char* description) {
    if (!std::isfinite(input)) {
        throw std::invalid_argument(std::string("Non-finite 3D thermal checkpoint ") + description);
    }
}

} // namespace

void write_transient_heat_checkpoint_3d(
    std::ostream& output,
    const fem::TetrahedralMesh3D& mesh,
    double time,
    const std::vector<double>& temperature
) {
    mesh.validate();
    check_finite(time, "time");
    if (temperature.size() != mesh.nodes.size()) {
        throw std::invalid_argument("3D thermal checkpoint temperature count must match mesh nodes.");
    }
    for (double temperature_value : temperature) {
        check_finite(temperature_value, "temperature");
    }
    if (!output.good()) {
        throw std::runtime_error("3D thermal checkpoint output stream is not writable.");
    }

    // Stage the complete serialization: invalid inputs never write a partial
    // checkpoint and the caller's output formatting remains untouched.
    std::ostringstream staged;
    staged.imbue(std::locale::classic());
    staged << std::setprecision(std::numeric_limits<double>::max_digits10);
    staged << "NOVA_THERMAL_CHECKPOINT_3D 1\n";
    staged << "time " << time << '\n';
    staged << "nodes " << mesh.nodes.size() << '\n';
    for (const auto& node : mesh.nodes) {
        staged << node.x << ' ' << node.y << ' ' << node.z << '\n';
    }
    staged << "elements " << mesh.elements.size() << '\n';
    for (const auto& element : mesh.elements) {
        staged << element[0] << ' ' << element[1] << ' '
               << element[2] << ' ' << element[3] << '\n';
    }
    staged << "temperature " << temperature.size() << '\n';
    for (double temperature_value : temperature) {
        staged << temperature_value << '\n';
    }
    staged << "end\n";
    if (!staged) {
        throw std::runtime_error("Failed to serialize 3D thermal checkpoint.");
    }
    output << staged.str();
    if (!output) {
        throw std::runtime_error("Failed to write 3D thermal checkpoint.");
    }
}

TransientHeatCheckpoint3D read_transient_heat_checkpoint_3d(
    std::istream& input,
    const fem::TetrahedralMesh3D& expected_mesh
) {
    expected_mesh.validate();
    expect(input, "NOVA_THERMAL_CHECKPOINT_3D");
    if (value<unsigned>(input, "version") != 1U) {
        throw std::invalid_argument("Unsupported 3D thermal checkpoint version.");
    }
    expect(input, "time");
    TransientHeatCheckpoint3D checkpoint;
    checkpoint.time = value<double>(input, "time");
    check_finite(checkpoint.time, "time");

    expect(input, "nodes");
    if (value<std::size_t>(input, "node count") != expected_mesh.nodes.size()) {
        throw std::invalid_argument("3D thermal checkpoint mesh node count mismatch.");
    }
    for (const auto& node : expected_mesh.nodes) {
        const double x = value<double>(input, "node x");
        const double y = value<double>(input, "node y");
        const double z = value<double>(input, "node z");
        check_finite(x, "node coordinate");
        check_finite(y, "node coordinate");
        check_finite(z, "node coordinate");
        if (x != node.x || y != node.y || z != node.z) {
            throw std::invalid_argument("3D thermal checkpoint mesh coordinates mismatch.");
        }
    }

    expect(input, "elements");
    if (value<std::size_t>(input, "element count") != expected_mesh.elements.size()) {
        throw std::invalid_argument("3D thermal checkpoint mesh element count mismatch.");
    }
    for (const auto& element : expected_mesh.elements) {
        for (std::size_t vertex = 0; vertex < element.size(); ++vertex) {
            if (value<std::size_t>(input, "element vertex") != element[vertex]) {
                throw std::invalid_argument("3D thermal checkpoint mesh connectivity mismatch.");
            }
        }
    }

    expect(input, "temperature");
    if (value<std::size_t>(input, "temperature count") != expected_mesh.nodes.size()) {
        throw std::invalid_argument("3D thermal checkpoint nodal temperature count mismatch.");
    }
    checkpoint.temperature.reserve(expected_mesh.nodes.size());
    for (std::size_t node = 0; node < expected_mesh.nodes.size(); ++node) {
        const double temperature = value<double>(input, "nodal temperature");
        check_finite(temperature, "temperature");
        checkpoint.temperature.push_back(temperature);
    }

    expect(input, "end");
    input >> std::ws;
    if (input.peek() != std::char_traits<char>::eof()) {
        throw std::invalid_argument("3D thermal checkpoint has trailing data.");
    }
    if (input.bad()) {
        throw std::runtime_error("Failed while reading 3D thermal checkpoint.");
    }
    return checkpoint;
}

} // namespace nova::io
