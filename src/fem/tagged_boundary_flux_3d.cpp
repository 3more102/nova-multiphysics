#include "nova/fem/tagged_boundary_flux_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace nova::fem {
namespace {

struct Owner {
    std::size_t cell{};
    std::size_t opposite{};
    unsigned count{};
};

struct Totals {
    std::size_t faces{};
    long double area{};
    long double rate{};
};

BoundaryFace3D face_key(BoundaryFace3D face) {
    std::sort(face.begin(), face.end());
    return face;
}

[[nodiscard]] double finite_double(long double value) {
    if (!std::isfinite(value) ||
        std::abs(value) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Boundary heat-flow diagnostic is outside finite double range.");
    }
    return static_cast<double>(value);
}

} // namespace

std::vector<TaggedBoundaryFlux3D> evaluate_tagged_boundary_flux_3d(
    const TetrahedralMesh3D& mesh,
    const std::span<const BoundaryFace3D> exterior_faces,
    const std::span<const int> physical_tags,
    const std::span<const P1TetrahedralFlux3D> cell_fields
) {
    mesh.validate();
    if (exterior_faces.size() != physical_tags.size()) {
        throw std::invalid_argument("Every boundary face needs exactly one physical tag.");
    }
    if (cell_fields.size() != mesh.elements.size()) {
        throw std::invalid_argument("Boundary flux requires one P1 flux per tetrahedron.");
    }
    for (const auto& field : cell_fields) {
        if (!std::isfinite(field.flux.x) || !std::isfinite(field.flux.y) ||
            !std::isfinite(field.flux.z)) {
            throw std::invalid_argument("Boundary flux vectors must be finite.");
        }
    }

    // A face shared by two cells is interior, even if the caller labels it.
    // Keep the opposite vertex to orient face normals independently of Gmsh
    // record winding and tetrahedral element ordering.
    std::map<BoundaryFace3D, Owner> owners;
    for (std::size_t cell = 0; cell < mesh.elements.size(); ++cell) {
        const auto& tetra = mesh.elements[cell];
        for (std::size_t omitted = 0; omitted < tetra.size(); ++omitted) {
            BoundaryFace3D face{};
            std::size_t at = 0;
            for (std::size_t j = 0; j < tetra.size(); ++j) {
                if (j != omitted) face[at++] = tetra[j];
            }
            auto [it, inserted] = owners.try_emplace(
                face_key(face), Owner{cell, tetra[omitted], 1U}
            );
            if (!inserted) ++it->second.count;
        }
    }

    std::set<BoundaryFace3D> visited;
    std::map<int, Totals> totals;
    for (std::size_t index = 0; index < exterior_faces.size(); ++index) {
        const int tag = physical_tags[index];
        if (tag < 0) throw std::invalid_argument("Boundary physical tags must be nonnegative.");
        const auto face = face_key(exterior_faces[index]);
        if (face[0] == face[1] || face[1] == face[2] ||
            face[2] >= mesh.nodes.size()) {
            throw std::invalid_argument("Invalid boundary triangle vertex indices.");
        }
        if (!visited.insert(face).second) {
            throw std::invalid_argument("Duplicate exterior face in boundary flux report.");
        }
        const auto found = owners.find(face);
        if (found == owners.end() || found->second.count != 1U) {
            throw std::invalid_argument("Boundary heat-flow face is not an exterior mesh face.");
        }

        const auto& a = mesh.nodes[face[0]];
        const auto& b = mesh.nodes[face[1]];
        const auto& c = mesh.nodes[face[2]];
        const auto& other = mesh.nodes[found->second.opposite];
        const std::array<long double, 3> ab{
            static_cast<long double>(b.x) - a.x,
            static_cast<long double>(b.y) - a.y,
            static_cast<long double>(b.z) - a.z
        };
        const std::array<long double, 3> ac{
            static_cast<long double>(c.x) - a.x,
            static_cast<long double>(c.y) - a.y,
            static_cast<long double>(c.z) - a.z
        };
        // Vector magnitude is the triangle's area; initial direction depends
        // on the sorted vertex winding, then is flipped away from the tetra.
        std::array<long double, 3> area_vector{
            (ab[1] * ac[2] - ab[2] * ac[1]) / 2,
            (ab[2] * ac[0] - ab[0] * ac[2]) / 2,
            (ab[0] * ac[1] - ab[1] * ac[0]) / 2
        };
        const long double inward =
            area_vector[0] * (static_cast<long double>(other.x) - a.x) +
            area_vector[1] * (static_cast<long double>(other.y) - a.y) +
            area_vector[2] * (static_cast<long double>(other.z) - a.z);
        if (!std::isfinite(inward) || inward == 0) {
            throw std::overflow_error("Boundary face orientation is not representable.");
        }
        if (inward > 0) {
            for (auto& component : area_vector) component = -component;
        }

        const long double area = std::hypot(
            area_vector[0], area_vector[1], area_vector[2]
        );
        const auto& flux = cell_fields[found->second.cell].flux;
        const long double rate =
            area_vector[0] * flux.x + area_vector[1] * flux.y +
            area_vector[2] * flux.z;
        if (!std::isfinite(area) || !std::isfinite(rate) || !(area > 0)) {
            throw std::overflow_error("Boundary face area or heat rate is non-finite.");
        }
        auto& group = totals[tag];
        ++group.faces;
        group.area += area;
        group.rate += rate;
    }

    std::vector<TaggedBoundaryFlux3D> result;
    result.reserve(totals.size());
    for (const auto& [tag, value] : totals) {
        result.push_back({tag, value.faces,
                          finite_double(value.area), finite_double(value.rate)});
    }
    return result;
}

} // namespace nova::fem
