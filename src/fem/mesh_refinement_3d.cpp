#include "nova/fem/mesh_refinement_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using EdgeKey = std::array<std::size_t, 2>;
using FaceKey = std::array<std::size_t, 3>;

[[nodiscard]] EdgeKey edge_key(const std::size_t a, const std::size_t b) {
    return {std::min(a, b), std::max(a, b)};
}

[[nodiscard]] FaceKey face_key(FaceKey face) {
    std::sort(face.begin(), face.end());
    return face;
}

[[nodiscard]] bool same_winding(const BoundaryFace3D& a, const BoundaryFace3D& b) {
    return a == b ||
        (a[0] == b[1] && a[1] == b[2]) ||
        (a[0] == b[2] && a[1] == b[0]);
}

[[nodiscard]] double squared_distance(const Point3D& a, const Point3D& b) {
    // hypot avoids intermediate squaring overflow for ordinary finite deltas;
    // numeric range beyond representable edge lengths is unsupported.
    const double length = std::hypot(a.x - b.x, a.y - b.y, a.z - b.z);
    return length * length;
}

} // namespace

RefinedTetrahedralMesh3D refine_tetrahedra_uniformly_3d(
    const TetrahedralMesh3D& mesh,
    const std::vector<BoundaryFace3D>& boundary_faces
) {
    mesh.validate();
    const auto exterior = extract_boundary_faces_3d(mesh);
    std::map<FaceKey, BoundaryFace3D> exterior_faces;
    for (const auto& face : exterior) {
        exterior_faces.emplace(face_key(face), face);
    }
    std::set<FaceKey> used_faces;
    for (const auto& face : boundary_faces) {
        for (const auto node : face) {
            if (node >= mesh.nodes.size()) {
                throw std::out_of_range("Refinement boundary face references a missing node.");
            }
        }
        const auto found = exterior_faces.find(face_key(face));
        if (found == exterior_faces.end() || !same_winding(face, found->second)) {
            throw std::invalid_argument(
                "Refinement boundary faces must be distinct outward-oriented exterior faces."
            );
        }
        if (!used_faces.insert(face_key(face)).second) {
            throw std::invalid_argument("Refinement boundary face is duplicated.");
        }
    }

    std::set<EdgeKey> edges;
    for (const auto& t : mesh.elements) {
        for (std::size_t i = 0; i < 4U; ++i) {
            for (std::size_t j = i + 1U; j < 4U; ++j) {
                edges.insert(edge_key(t[i], t[j]));
            }
        }
    }

    RefinedTetrahedralMesh3D result;
    result.coarse_node_count = mesh.nodes.size();
    result.coarse_element_count = mesh.elements.size();
    if (edges.size() > result.mesh.nodes.max_size() - mesh.nodes.size() ||
        edges.size() > result.midpoint_parent_edges.max_size() ||
        mesh.elements.size() > result.mesh.elements.max_size() / 8U ||
        mesh.elements.size() > result.parent_elements.max_size() / 8U ||
        boundary_faces.size() > result.boundary_faces.max_size() / 4U ||
        boundary_faces.size() > result.parent_boundary_faces.max_size() / 4U) {
        throw std::length_error("Uniform tetrahedral refinement exceeds capacity.");
    }

    result.mesh.nodes = mesh.nodes;
    result.mesh.nodes.reserve(mesh.nodes.size() + edges.size());
    result.midpoint_parent_edges.reserve(edges.size());
    result.mesh.elements.reserve(8U * mesh.elements.size());
    result.parent_elements.reserve(8U * mesh.elements.size());
    result.boundary_faces.reserve(4U * boundary_faces.size());
    result.parent_boundary_faces.reserve(4U * boundary_faces.size());

    // Ordered canonical edges make generated node indices independent of
    // tetrahedron insertion order and local element permutation.
    std::map<EdgeKey, std::size_t> midpoint_ids;
    for (const auto& edge : edges) {
        const auto& a = mesh.nodes[edge[0]];
        const auto& b = mesh.nodes[edge[1]];
        const Point3D midpoint{
            0.5 * a.x + 0.5 * b.x,
            0.5 * a.y + 0.5 * b.y,
            0.5 * a.z + 0.5 * b.z
        };
        if (!std::isfinite(midpoint.x) || !std::isfinite(midpoint.y) ||
            !std::isfinite(midpoint.z)) {
            throw std::overflow_error("Refined tetrahedron midpoint is non-finite.");
        }
        midpoint_ids.emplace(edge, result.mesh.nodes.size());
        result.mesh.nodes.push_back(midpoint);
        result.midpoint_parent_edges.push_back(edge);
    }

    const auto midpoint_id = [&midpoint_ids](const std::size_t a, const std::size_t b) {
        return midpoint_ids.at(edge_key(a, b));
    };
    for (std::size_t parent = 0; parent < mesh.elements.size(); ++parent) {
        const auto& t = mesh.elements[parent];
        const std::size_t a = t[0], b = t[1], c = t[2], d = t[3];
        const auto ab = midpoint_id(a, b), ac = midpoint_id(a, c);
        const auto ad = midpoint_id(a, d), bc = midpoint_id(b, c);
        const auto bd = midpoint_id(b, d), cd = midpoint_id(c, d);
        const std::array<Tetrahedron3D, 4> corners{{
            {{a, ab, ac, ad}}, {{b, ab, bc, bd}},
            {{c, ac, bc, cd}}, {{d, ad, bd, cd}}
        }};
        for (const auto& tet : corners) {
            result.mesh.elements.push_back(tet);
            result.parent_elements.push_back(parent);
        }

        // Split the central octahedron into four tetrahedra around one
        // interior diagonal. The six midpoint nodes lie on the octahedron.
        // Its boundary is unchanged for all three possible choices.
        const std::array<std::array<std::size_t, 6>, 3> diagonals{{
            {{ab, cd, ac, ad, bd, bc}},
            {{ac, bd, ab, ad, cd, bc}},
            {{ad, bc, ab, ac, cd, bd}}
        }};
        std::size_t best = 0;
        double best_squared_length = std::numeric_limits<double>::infinity();
        for (std::size_t choice = 0; choice < diagonals.size(); ++choice) {
            const auto& diagonal = diagonals[choice];
            const double squared_length = squared_distance(
                result.mesh.nodes[diagonal[0]], result.mesh.nodes[diagonal[1]]
            );
            if (!std::isfinite(squared_length)) {
                throw std::overflow_error("Central octahedron diagonal length overflows.");
            }
            if (squared_length < best_squared_length) {
                best_squared_length = squared_length;
                best = choice;
            }
        }
        const auto& v = diagonals[best];
        // v[2..5] winds around the chosen diagonal in its octahedron.
        for (std::size_t i = 0; i < 4U; ++i) {
            result.mesh.elements.push_back({
                v[0], v[1], v[2U + i], v[2U + (i + 1U) % 4U]
            });
            result.parent_elements.push_back(parent);
        }
    }

    for (std::size_t parent = 0; parent < boundary_faces.size(); ++parent) {
        const auto& f = boundary_faces[parent];
        const auto a = f[0], b = f[1], c = f[2];
        const auto ab = midpoint_id(a, b);
        const auto bc = midpoint_id(b, c);
        const auto ca = midpoint_id(c, a);
        const std::array<BoundaryFace3D, 4> children{{
            {{a, ab, ca}}, {{ab, b, bc}},
            {{ca, bc, c}}, {{ab, bc, ca}}
        }};
        for (const auto& child : children) {
            result.boundary_faces.push_back(child);
            result.parent_boundary_faces.push_back(parent);
        }
    }

    result.mesh.validate();
    // This checks that the user-supplied face subdivision exactly matches
    // the refined mesh's exterior rather than merely being geometrically close.
    if (!boundary_faces.empty()) {
        std::set<FaceKey> refined_exterior;
        for (const auto& face : extract_boundary_faces_3d(result.mesh)) {
            refined_exterior.insert(face_key(face));
        }
        for (const auto& face : result.boundary_faces) {
            if (!refined_exterior.contains(face_key(face))) {
                throw std::runtime_error("Refined boundary face is not exterior.");
            }
        }
    }
    return result;
}

std::vector<double> prolongate_p1_uniform_3d(
    const RefinedTetrahedralMesh3D& refined,
    const std::span<const double> coarse_nodal_values
) {
    if (coarse_nodal_values.size() != refined.coarse_node_count ||
        refined.coarse_node_count > refined.mesh.nodes.size() ||
        refined.mesh.nodes.size() - refined.coarse_node_count !=
            refined.midpoint_parent_edges.size()) {
        throw std::invalid_argument("3D P1 prolongation has inconsistent DOF sizes.");
    }
    std::vector<double> result(coarse_nodal_values.begin(), coarse_nodal_values.end());
    for (const double value : result) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Coarse P1 values must be finite.");
        }
    }
    result.reserve(refined.mesh.nodes.size());
    for (const auto& edge : refined.midpoint_parent_edges) {
        if (edge[0] >= refined.coarse_node_count ||
            edge[1] >= refined.coarse_node_count) {
            throw std::out_of_range("3D P1 prolongation references missing parent DOF.");
        }
        const double interpolated = 0.5 * result[edge[0]] + 0.5 * result[edge[1]];
        if (!std::isfinite(interpolated)) {
            throw std::overflow_error("3D P1 interpolated value is non-finite.");
        }
        result.push_back(interpolated);
    }
    return result;
}

std::vector<double> prolongate_cellwise_constant_uniform_3d(
    const RefinedTetrahedralMesh3D& refined,
    const std::span<const double> coarse_cell_values
) {
    if (coarse_cell_values.size() != refined.coarse_element_count ||
        refined.mesh.elements.size() != refined.parent_elements.size()) {
        throw std::invalid_argument("3D cell transfer has inconsistent dimensions.");
    }
    for (const double value : coarse_cell_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Coarse cell values must be finite.");
        }
    }
    std::vector<double> result;
    result.reserve(refined.parent_elements.size());
    for (const auto parent : refined.parent_elements) {
        if (parent >= coarse_cell_values.size()) {
            throw std::out_of_range("Cell transfer references missing parent element.");
        }
        result.push_back(coarse_cell_values[parent]);
    }
    return result;
}

} // namespace nova::fem
