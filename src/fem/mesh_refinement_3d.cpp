#include "nova/fem/mesh_refinement_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace nova::fem {
namespace {
using Edge = std::array<std::size_t, 2>;

[[nodiscard]] Edge canonical_edge(std::size_t a, std::size_t b) {
    if (a > b) std::swap(a, b);
    return {a, b};
}

[[nodiscard]] Point3D midpoint(const Point3D& a, const Point3D& b) {
    // Half differences avoid overflow in a+b for large, same-signed inputs.
    return {
        a.x + (b.x - a.x) * 0.5,
        a.y + (b.y - a.y) * 0.5,
        a.z + (b.z - a.z) * 0.5
    };
}

} // namespace

UniformTetrahedralRefinement3D refine_tetrahedra_uniform_3d(
    const TetrahedralMesh3D& coarse
) {
    coarse.validate();
    UniformTetrahedralRefinement3D result;
    result.original_node_count = coarse.nodes.size();
    result.mesh.nodes = coarse.nodes;
    if (coarse.elements.size() > result.mesh.elements.max_size() / 8U ||
        coarse.elements.size() > result.parent_elements.max_size() / 8U) {
        throw std::length_error("3D refinement exceeds tetrahedron storage capacity.");
    }
    result.mesh.elements.reserve(coarse.elements.size() * 8U);
    result.parent_elements.reserve(coarse.elements.size() * 8U);
    std::map<Edge, std::size_t> shared_midpoints;

    auto get_midpoint = [&](std::size_t a, std::size_t b) {
        const Edge edge = canonical_edge(a, b);
        if (const auto it = shared_midpoints.find(edge); it != shared_midpoints.end()) {
            return it->second;
        }
        if (result.mesh.nodes.size() == result.mesh.nodes.max_size()) {
            throw std::length_error("3D refinement exceeds vertex storage capacity.");
        }
        const auto new_node = result.mesh.nodes.size();
        const auto value = midpoint(coarse.nodes[a], coarse.nodes[b]);
        if (!std::isfinite(value.x) || !std::isfinite(value.y) ||
            !std::isfinite(value.z)) {
            throw std::overflow_error("3D edge midpoint is nonfinite.");
        }
        result.mesh.nodes.push_back(value);
        result.midpoint_parent_edges.push_back(edge);
        shared_midpoints.emplace(edge, new_node);
        return new_node;
    };

    for (std::size_t parent = 0; parent < coarse.elements.size(); ++parent) {
        const auto& t = coarse.elements[parent];
        const auto ab = get_midpoint(t[0], t[1]);
        const auto ac = get_midpoint(t[0], t[2]);
        const auto ad = get_midpoint(t[0], t[3]);
        const auto bc = get_midpoint(t[1], t[2]);
        const auto bd = get_midpoint(t[1], t[3]);
        const auto cd = get_midpoint(t[2], t[3]);
        const std::array<Tetrahedron3D, 8> children{{
            {{t[0], ab, ac, ad}},
            {{t[1], ab, bc, bd}},
            {{t[2], ac, bc, cd}},
            {{t[3], ad, bd, cd}},
            {{ab, ac, ad, cd}},
            {{ab, ad, bd, cd}},
            {{ab, bd, bc, cd}},
            {{ab, bc, ac, cd}}
        }};
        for (const auto& child : children) {
            result.mesh.elements.push_back(child);
            result.parent_elements.push_back(parent);
        }
    }
    result.mesh.validate();

    const auto exterior = extract_boundary_faces_3d(coarse);
    if (exterior.size() > result.boundary_faces.max_size() / 4U ||
        exterior.size() > result.parent_boundary_faces.max_size() / 4U) {
        throw std::length_error("3D refinement exceeds boundary-face storage capacity.");
    }
    result.boundary_faces.reserve(exterior.size() * 4U);
    result.parent_boundary_faces.reserve(exterior.size() * 4U);
    for (std::size_t parent = 0; parent < exterior.size(); ++parent) {
        const auto& f = exterior[parent];
        const auto ab = shared_midpoints.at(canonical_edge(f[0], f[1]));
        const auto bc = shared_midpoints.at(canonical_edge(f[1], f[2]));
        const auto ca = shared_midpoints.at(canonical_edge(f[2], f[0]));
        // Retain parent exterior winding, hence outward physical orientation.
        const std::array<BoundaryFace3D, 4> children{{
            {{f[0], ab, ca}}, {{ab, f[1], bc}},
            {{ca, bc, f[2]}}, {{ab, bc, ca}}
        }};
        for (const auto& face : children) {
            result.boundary_faces.push_back(face);
            result.parent_boundary_faces.push_back(parent);
        }
    }
    return result;
}

std::vector<double> prolongate_tetrahedral_p1_3d(
    const UniformTetrahedralRefinement3D& refinement,
    const std::vector<double>& coarse_values
) {
    if (coarse_values.size() != refinement.original_node_count ||
        refinement.original_node_count > refinement.mesh.nodes.size() ||
        refinement.midpoint_parent_edges.size() !=
            refinement.mesh.nodes.size() - refinement.original_node_count) {
        throw std::invalid_argument("3D P1 transfer field or refinement dimensions mismatch.");
    }
    for (const auto value : coarse_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("3D P1 transfer requires finite nodal values.");
        }
    }
    std::vector<double> fine = coarse_values;
    fine.reserve(refinement.mesh.nodes.size());
    for (const Edge& endpoints : refinement.midpoint_parent_edges) {
        if (endpoints[0] >= refinement.original_node_count ||
            endpoints[1] >= refinement.original_node_count ||
            endpoints[0] == endpoints[1]) {
            throw std::invalid_argument("3D P1 transfer has invalid edge metadata.");
        }
        const double value = 0.5 * coarse_values[endpoints[0]] +
                             0.5 * coarse_values[endpoints[1]];
        if (!std::isfinite(value)) {
            throw std::overflow_error("3D P1 interpolated midpoint is nonfinite.");
        }
        fine.push_back(value);
    }
    return fine;
}

} // namespace nova::fem
