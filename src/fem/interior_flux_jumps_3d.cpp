#include "nova/fem/interior_flux_jumps_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace nova::fem {
namespace {

struct FaceOwners {
    std::size_t first_cell{};
    std::size_t second_cell{};
    std::size_t first_opposite{};
    unsigned occurrences{};
};

struct Vector3L {
    long double x{};
    long double y{};
    long double z{};
};

[[nodiscard]] Vector3L difference(const Point3D& a, const Point3D& b) {
    return {
        static_cast<long double>(a.x) - b.x,
        static_cast<long double>(a.y) - b.y,
        static_cast<long double>(a.z) - b.z
    };
}

[[nodiscard]] Vector3L cross(const Vector3L& a, const Vector3L& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

[[nodiscard]] long double dot(const Vector3L& a, const Vector3L& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

[[nodiscard]] long double length(const Vector3L& value) {
    return std::hypot(value.x, value.y, value.z);
}

[[nodiscard]] double checked_double(long double value) {
    if (!std::isfinite(value) ||
        std::abs(value) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Interior flux-jump diagnostic exceeds the finite range.");
    }
    return static_cast<double>(value);
}

} // namespace

InteriorFluxJumpReport3D evaluate_interior_flux_jumps_3d(
    const TetrahedralMesh3D& mesh,
    const std::span<const P1TetrahedralFlux3D> cell_fields
) {
    mesh.validate();
    if (cell_fields.size() != mesh.elements.size()) {
        throw std::invalid_argument("3D jump diagnostics require one flux per tetrahedron.");
    }
    for (const auto& field : cell_fields) {
        if (!std::isfinite(field.flux.x) || !std::isfinite(field.flux.y) ||
            !std::isfinite(field.flux.z)) {
            throw std::invalid_argument("3D jump diagnostics require finite cell fluxes.");
        }
    }

    std::map<BoundaryFace3D, FaceOwners> faces;
    for (std::size_t cell = 0; cell < mesh.elements.size(); ++cell) {
        const auto& tetra = mesh.elements[cell];
        for (std::size_t omitted = 0; omitted < 4; ++omitted) {
            BoundaryFace3D key{};
            std::size_t index = 0;
            for (std::size_t vertex = 0; vertex < 4; ++vertex) {
                if (vertex != omitted) key[index++] = tetra[vertex];
            }
            std::sort(key.begin(), key.end());
            const auto [at, inserted] = faces.try_emplace(
                key, FaceOwners{cell, 0, tetra[omitted], 1U}
            );
            if (!inserted) {
                // mesh.validate() already rejected nonmanifold faces.
                at->second.second_cell = cell;
                ++at->second.occurrences;
            }
        }
    }

    InteriorFluxJumpReport3D result;
    result.cell_indicators.resize(mesh.elements.size(), 0.0);
    std::vector<long double> squared_cell(mesh.elements.size(), 0.0L);
    long double squared_global = 0.0L;
    long double maximum_jump = 0.0L;
    for (const auto& [face, owners] : faces) {
        if (owners.occurrences == 1U) continue;
        const auto& a = mesh.nodes[face[0]];
        const auto& b = mesh.nodes[face[1]];
        const auto& c = mesh.nodes[face[2]];
        const auto& opposite = mesh.nodes[owners.first_opposite];
        const auto ab = difference(b, a);
        const auto ac = difference(c, a);
        auto area_vector = cross(ab, ac);
        area_vector.x *= 0.5L;
        area_vector.y *= 0.5L;
        area_vector.z *= 0.5L;
        const long double inward = dot(area_vector, difference(opposite, a));
        if (!std::isfinite(inward) || inward == 0.0L) {
            throw std::overflow_error("Interior face orientation cannot be represented.");
        }
        if (inward > 0.0L) {
            area_vector.x = -area_vector.x;
            area_vector.y = -area_vector.y;
            area_vector.z = -area_vector.z;
        }
        const long double area = length(area_vector);
        const long double diameter = std::max({
            length(ab), length(ac), length(difference(c, b))
        });
        const auto& first = cell_fields[owners.first_cell].flux;
        const auto& second = cell_fields[owners.second_cell].flux;
        const Vector3L flux_difference{
            static_cast<long double>(first.x) - second.x,
            static_cast<long double>(first.y) - second.y,
            static_cast<long double>(first.z) - second.z
        };
        const long double jump = dot(flux_difference, area_vector) / area;
        const long double contribution = diameter * area * jump * jump;
        if (!(area > 0.0L) || !(diameter > 0.0L) ||
            !std::isfinite(jump) || !std::isfinite(contribution)) {
            throw std::overflow_error("Interior flux jump is not representable.");
        }
        squared_global += contribution;
        squared_cell[owners.first_cell] += 0.5L * contribution;
        squared_cell[owners.second_cell] += 0.5L * contribution;
        maximum_jump = std::max(maximum_jump, std::abs(jump));
        result.faces.push_back({
            face, owners.first_cell, owners.second_cell,
            checked_double(area), checked_double(diameter),
            checked_double(jump)
        });
    }
    result.global_indicator = checked_double(std::sqrt(squared_global));
    result.maximum_absolute_normal_jump = checked_double(maximum_jump);
    for (std::size_t i = 0; i < squared_cell.size(); ++i) {
        result.cell_indicators[i] = checked_double(std::sqrt(squared_cell[i]));
    }
    return result;
}

} // namespace nova::fem
