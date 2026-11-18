#include "nova/fem/mesh3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace nova::fem {
namespace {

struct Vector3 {
    double x{};
    double y{};
    double z{};
};

using FaceKey = std::array<std::size_t, 3>;

[[nodiscard]] Vector3 subtract(const Point3D& a, const Point3D& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

[[nodiscard]] Vector3 cross(const Vector3& a, const Vector3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

[[nodiscard]] double dot(const Vector3& a, const Vector3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

[[nodiscard]] double norm(const Vector3& value) {
    return std::hypot(value.x, value.y, value.z);
}

[[nodiscard]] Vector3 divide(const Vector3& value, double scale) {
    return {value.x / scale, value.y / scale, value.z / scale};
}

[[nodiscard]] double determinant(
    const Point3D& p0,
    const Point3D& p1,
    const Point3D& p2,
    const Point3D& p3
) {
    const Vector3 e1 = subtract(p1, p0);
    const Vector3 e2 = subtract(p2, p0);
    const Vector3 e3 = subtract(p3, p0);
    return dot(e1, cross(e2, e3));
}

[[nodiscard]] double edge_length(
    const Point3D& a,
    const Point3D& b
) {
    return norm(subtract(a, b));
}

[[nodiscard]] FaceKey canonical_face(BoundaryFace3D face) {
    std::sort(face.begin(), face.end());
    return face;
}

struct OrientedFace {
    BoundaryFace3D nodes{};
    std::size_t opposite{};
};

[[nodiscard]] std::array<OrientedFace, 4> outward_faces(
    const TetrahedralMesh3D& mesh,
    const Tetrahedron3D& tetrahedron
) {
    std::array<OrientedFace, 4> faces{{
        {{{tetrahedron[1], tetrahedron[2], tetrahedron[3]}}, tetrahedron[0]},
        {{{tetrahedron[0], tetrahedron[3], tetrahedron[2]}}, tetrahedron[1]},
        {{{tetrahedron[0], tetrahedron[1], tetrahedron[3]}}, tetrahedron[2]},
        {{{tetrahedron[0], tetrahedron[2], tetrahedron[1]}}, tetrahedron[3]}
    }};

    for (OrientedFace& face : faces) {
        const Point3D& p0 = mesh.nodes[face.nodes[0]];
        const Point3D& p1 = mesh.nodes[face.nodes[1]];
        const Point3D& p2 = mesh.nodes[face.nodes[2]];
        const Point3D& opposite = mesh.nodes[face.opposite];
        const Vector3 normal = cross(
            subtract(p1, p0), subtract(p2, p0)
        );
        if (dot(normal, subtract(opposite, p0)) > 0.0) {
            std::swap(face.nodes[1], face.nodes[2]);
        }
    }
    return faces;
}

[[nodiscard]] std::array<double, 6> tetrahedron_edge_lengths(
    const TetrahedralMesh3D& mesh,
    const Tetrahedron3D& tetrahedron
) {
    return {{
        edge_length(mesh.nodes[tetrahedron[0]], mesh.nodes[tetrahedron[1]]),
        edge_length(mesh.nodes[tetrahedron[0]], mesh.nodes[tetrahedron[2]]),
        edge_length(mesh.nodes[tetrahedron[0]], mesh.nodes[tetrahedron[3]]),
        edge_length(mesh.nodes[tetrahedron[1]], mesh.nodes[tetrahedron[2]]),
        edge_length(mesh.nodes[tetrahedron[1]], mesh.nodes[tetrahedron[3]]),
        edge_length(mesh.nodes[tetrahedron[2]], mesh.nodes[tetrahedron[3]])
    }};
}

[[nodiscard]] TetrahedronQuality3D evaluate_tetrahedron_quality(
    const TetrahedralMesh3D& mesh,
    const Tetrahedron3D& tetrahedron
) {
    const Point3D& p0 = mesh.nodes[tetrahedron[0]];
    const Point3D& p1 = mesh.nodes[tetrahedron[1]];
    const Point3D& p2 = mesh.nodes[tetrahedron[2]];
    const Point3D& p3 = mesh.nodes[tetrahedron[3]];
    const double det = determinant(p0, p1, p2, p3);
    const double volume = std::abs(det) / 6.0;
    const auto lengths = tetrahedron_edge_lengths(mesh, tetrahedron);

    const double maximum_edge = *std::max_element(lengths.begin(), lengths.end());
    double edge_square_sum = 0.0;
    for (const double length : lengths) {
        const double relative_length = length / maximum_edge;
        edge_square_sum += relative_length * relative_length;
    }
    const auto e1 = divide(subtract(p1, p0), maximum_edge);
    const auto e2 = divide(subtract(p2, p0), maximum_edge);
    const auto e3 = divide(subtract(p3, p0), maximum_edge);
    const double scaled_det = std::abs(dot(e1, cross(e2, e3)));
    const double mean_ratio =
        12.0 * std::pow(scaled_det / 2.0, 2.0 / 3.0) / edge_square_sum;

    const auto faces = outward_faces(mesh, tetrahedron);
    std::array<Vector3, 4> unit_normals{};
    for (std::size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i].nodes;
        const auto& origin = mesh.nodes[face[0]];
        const auto a = divide(subtract(mesh.nodes[face[1]], origin), maximum_edge);
        const auto b = divide(subtract(mesh.nodes[face[2]], origin), maximum_edge);
        const auto normal = cross(a, b);
        unit_normals[i] = divide(normal, norm(normal));
    }
    double minimum_dihedral = std::acos(-1.0);
    for (std::size_t i = 0; i < faces.size(); ++i) {
        for (std::size_t j = i + 1; j < faces.size(); ++j) {
            const double cosine = std::clamp(-dot(unit_normals[i], unit_normals[j]), -1.0, 1.0);
            minimum_dihedral = std::min(minimum_dihedral, std::acos(cosine));
        }
    }

    return {
        det / 6.0,
        volume,
        mean_ratio,
        minimum_dihedral,
        maximum_edge
    };
}

struct FaceRecord {
    BoundaryFace3D outward{};
    bool interior{};
};

// Both faces have the same three vertex IDs. Cyclic permutations preserve winding.
[[nodiscard]] bool same_winding(const BoundaryFace3D& a, const BoundaryFace3D& b) {
    return a == b || (a[0] == b[1] && a[1] == b[2]) ||
                     (a[0] == b[2] && a[1] == b[0]);
}

[[nodiscard]] std::map<FaceKey, FaceRecord> validated_faces(const TetrahedralMesh3D& mesh) {
    const auto& nodes = mesh.nodes;
    const auto& elements = mesh.elements;
    if (nodes.empty()) {
        throw std::invalid_argument("A 3D tetrahedral mesh must contain nodes.");
    }
    if (elements.empty()) {
        throw std::invalid_argument("A 3D tetrahedral mesh must contain elements.");
    }
    for (const Point3D& node : nodes) {
        if (!std::isfinite(node.x) || !std::isfinite(node.y) ||
            !std::isfinite(node.z)) {
            throw std::invalid_argument("3D mesh node coordinates must be finite.");
        }
    }

    constexpr double degeneracy_factor =
        128.0 * std::numeric_limits<double>::epsilon();
    std::set<std::array<std::size_t, 4>> unique_elements;
    std::map<FaceKey, FaceRecord> records;
    std::vector<bool> used(nodes.size(), false);

    for (const Tetrahedron3D& tetrahedron : elements) {
        for (const std::size_t node : tetrahedron) {
            if (node >= nodes.size()) {
                throw std::out_of_range(
                    "Tetrahedron connectivity references a missing node."
                );
            }
            used[node] = true;
        }
        for (std::size_t i = 0; i < 4U; ++i) {
            for (std::size_t j = i + 1U; j < 4U; ++j) {
                if (tetrahedron[i] == tetrahedron[j]) {
                    throw std::invalid_argument(
                        "A tetrahedron cannot repeat a node index."
                    );
                }
            }
        }

        auto canonical_element = tetrahedron;
        std::sort(canonical_element.begin(), canonical_element.end());
        if (!unique_elements.insert(canonical_element).second) {
            throw std::invalid_argument("Mesh contains a duplicate tetrahedron.");
        }

        const auto lengths = tetrahedron_edge_lengths(mesh, tetrahedron);
        double maximum_edge = 0.0;
        for (const double length : lengths) {
            maximum_edge = std::max(maximum_edge, length);
        }
        const double maximum_edge_cubed =
            maximum_edge * maximum_edge * maximum_edge;
        const double det = determinant(
            nodes[tetrahedron[0]], nodes[tetrahedron[1]],
            nodes[tetrahedron[2]], nodes[tetrahedron[3]]
        );
        if (!std::isfinite(maximum_edge_cubed) || !std::isfinite(det)) {
            throw std::invalid_argument(
                "Tetrahedron geometry exceeds the supported numeric range."
            );
        }
        if (std::abs(det) / 6.0 == 0.0 ||
            std::abs(det) <= degeneracy_factor * maximum_edge_cubed) {
            throw std::invalid_argument(
                "Mesh contains a degenerate or near-degenerate tetrahedron."
            );
        }

        for (const OrientedFace& face : outward_faces(mesh, tetrahedron)) {
            const auto [it, inserted] = records.try_emplace(
                canonical_face(face.nodes), FaceRecord{face.nodes, false}
            );
            if (!inserted) {
                if (it->second.interior) {
                    throw std::invalid_argument("More than two tetrahedra share a face.");
                }
                if (same_winding(it->second.outward, face.nodes)) {
                    throw std::invalid_argument("Tetrahedra occupy the same side of a shared face.");
                }
                it->second.interior = true;
            }
        }
    }
    if (std::find(used.begin(), used.end(), false) != used.end()) {
        throw std::invalid_argument("Every 3D mesh vertex must belong to a tetrahedron.");
    }
    return records;
}

} // namespace

void TetrahedralMesh3D::validate() const {
    static_cast<void>(validated_faces(*this));
}

std::vector<BoundaryFace3D> extract_boundary_faces_3d(
    const TetrahedralMesh3D& mesh
) {
    const auto records = validated_faces(mesh);
    std::vector<BoundaryFace3D> boundary_faces;
    for (const auto& [key, record] : records) {
        static_cast<void>(key);
        if (!record.interior) {
            boundary_faces.push_back(record.outward);
        }
    }
    return boundary_faces;
}

MeshQuality3D evaluate_mesh_quality_3d(const TetrahedralMesh3D& mesh) {
    mesh.validate();
    MeshQuality3D result;
    result.elements.reserve(mesh.elements.size());
    result.minimum_mean_ratio = 1.0;
    result.minimum_dihedral_angle = std::acos(-1.0);

    for (const Tetrahedron3D& tetrahedron : mesh.elements) {
        const TetrahedronQuality3D quality =
            evaluate_tetrahedron_quality(mesh, tetrahedron);
        result.elements.push_back(quality);
        result.total_volume += quality.volume;
        result.minimum_mean_ratio = std::min(
            result.minimum_mean_ratio, quality.mean_ratio
        );
        result.minimum_dihedral_angle = std::min(
            result.minimum_dihedral_angle, quality.minimum_dihedral_angle
        );
    }
    if (!std::isfinite(result.total_volume)) {
        throw std::overflow_error("Total tetrahedral mesh volume overflowed.");
    }
    return result;
}

StructuredBoxMesh3D make_box_mesh_3d(
    const double x_min, const double x_max,
    const double y_min, const double y_max,
    const double z_min, const double z_max,
    const std::size_t elements_x,
    const std::size_t elements_y,
    const std::size_t elements_z
) {
    if (!std::isfinite(x_min) || !std::isfinite(x_max) ||
        !std::isfinite(y_min) || !std::isfinite(y_max) ||
        !std::isfinite(z_min) || !std::isfinite(z_max) ||
        x_max <= x_min || y_max <= y_min || z_max <= z_min) {
        throw std::invalid_argument(
            "Box bounds must be finite and have positive width, height, and depth."
        );
    }
    if (elements_x == 0U || elements_y == 0U || elements_z == 0U) {
        throw std::invalid_argument("Box mesh subdivisions must be positive.");
    }

    const std::size_t maximum = std::numeric_limits<std::size_t>::max();
    if (elements_x == maximum || elements_y == maximum ||
        elements_z == maximum) {
        throw std::length_error("Box mesh dimensions exceed the supported index range.");
    }
    const std::size_t nodes_x = elements_x + 1U;
    const std::size_t nodes_y = elements_y + 1U;
    const std::size_t nodes_z = elements_z + 1U;
    if (nodes_x > maximum / nodes_y ||
        nodes_x * nodes_y > maximum / nodes_z ||
        elements_x > maximum / elements_y ||
        elements_x * elements_y > maximum / elements_z) {
        throw std::length_error("Box mesh size overflows the supported index range.");
    }
    const std::size_t node_count = nodes_x * nodes_y * nodes_z;
    const std::size_t cell_count = elements_x * elements_y * elements_z;
    if (cell_count > maximum / 6U) {
        throw std::length_error("Box tetrahedron count overflows the supported index range.");
    }

    StructuredBoxMesh3D result;
    result.mesh.nodes.reserve(node_count);
    result.mesh.elements.reserve(6U * cell_count);

    const double step_x = (x_max - x_min) / static_cast<double>(elements_x);
    const double step_y = (y_max - y_min) / static_cast<double>(elements_y);
    const double step_z = (z_max - z_min) / static_cast<double>(elements_z);
    for (std::size_t k = 0; k < nodes_z; ++k) {
        for (std::size_t j = 0; j < nodes_y; ++j) {
            for (std::size_t i = 0; i < nodes_x; ++i) {
                const double x = i == elements_x
                    ? x_max : x_min + static_cast<double>(i) * step_x;
                const double y = j == elements_y
                    ? y_max : y_min + static_cast<double>(j) * step_y;
                const double z = k == elements_z
                    ? z_max : z_min + static_cast<double>(k) * step_z;
                result.mesh.nodes.push_back({x, y, z});
            }
        }
    }

    const auto node_id = [nodes_x, nodes_y](
        const std::size_t i,
        const std::size_t j,
        const std::size_t k
    ) {
        return i + nodes_x * (j + nodes_y * k);
    };

    for (std::size_t k = 0; k < elements_z; ++k) {
        for (std::size_t j = 0; j < elements_y; ++j) {
            for (std::size_t i = 0; i < elements_x; ++i) {
                const std::size_t v000 = node_id(i, j, k);
                const std::size_t v100 = node_id(i + 1U, j, k);
                const std::size_t v010 = node_id(i, j + 1U, k);
                const std::size_t v110 = node_id(i + 1U, j + 1U, k);
                const std::size_t v001 = node_id(i, j, k + 1U);
                const std::size_t v101 = node_id(i + 1U, j, k + 1U);
                const std::size_t v011 = node_id(i, j + 1U, k + 1U);
                const std::size_t v111 = node_id(i + 1U, j + 1U, k + 1U);

                const std::array<Tetrahedron3D, 6> local{{
                    {{v000, v100, v110, v111}},
                    {{v000, v110, v010, v111}},
                    {{v000, v010, v011, v111}},
                    {{v000, v011, v001, v111}},
                    {{v000, v001, v101, v111}},
                    {{v000, v101, v100, v111}}
                }};
                for (Tetrahedron3D tetrahedron : local) {
                    const double det = determinant(
                        result.mesh.nodes[tetrahedron[0]],
                        result.mesh.nodes[tetrahedron[1]],
                        result.mesh.nodes[tetrahedron[2]],
                        result.mesh.nodes[tetrahedron[3]]
                    );
                    if (det < 0.0) {
                        std::swap(tetrahedron[2], tetrahedron[3]);
                    }
                    result.mesh.elements.push_back(tetrahedron);
                }
            }
        }
    }

    result.boundary_faces = extract_boundary_faces_3d(result.mesh);
    return result;
}

} // namespace nova::fem
