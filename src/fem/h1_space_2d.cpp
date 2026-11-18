#include "nova/fem/h1_space_2d.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>

namespace nova::fem {
namespace {

[[nodiscard]] EdgeKey2D canonical_edge(
    const std::size_t first,
    const std::size_t second
) {
    return std::minmax(first, second);
}

void validate_edge_nodes(
    const TriangleMesh2D& mesh,
    const std::size_t first,
    const std::size_t second
) {
    if (first >= mesh.nodes.size() || second >= mesh.nodes.size()) {
        throw std::out_of_range("H1 edge references a missing mesh node.");
    }
    if (first == second) {
        throw std::invalid_argument("An H1 edge must have two distinct nodes.");
    }
}

} // namespace

H1DofManager2D::H1DofManager2D(
    const TriangleMesh2D& mesh,
    const H1TriangleOrder order
)
    : order_(order),
      vertex_count_(mesh.nodes.size()) {
    mesh.validate();

    if (order != H1TriangleOrder::P1 && order != H1TriangleOrder::P2) {
        throw std::invalid_argument("Unsupported H1 triangle polynomial order.");
    }

    std::map<EdgeKey2D, std::size_t> edge_indices;
    for (const Triangle2D& element : mesh.elements) {
        const std::array<EdgeKey2D, 3> local_edges{{
            canonical_edge(element[0], element[1]),
            canonical_edge(element[1], element[2]),
            canonical_edge(element[2], element[0])
        }};

        for (const EdgeKey2D& edge : local_edges) {
            if (edge_indices.find(edge) == edge_indices.end()) {
                const std::size_t index = edge_keys_.size();
                edge_indices.emplace(edge, index);
                edge_keys_.push_back(edge);
            }
        }
    }

    dof_count_ = vertex_count_ +
        (order_ == H1TriangleOrder::P2 ? edge_keys_.size() : 0U);

    element_dofs_.reserve(mesh.elements.size());

    for (const Triangle2D& element : mesh.elements) {
        H1ElementDofMap2D mapping;
        mapping.count = triangle_h1_dof_count(order_);

        mapping.dofs[0] = element[0];
        mapping.dofs[1] = element[1];
        mapping.dofs[2] = element[2];

        if (order_ == H1TriangleOrder::P2) {
            const std::array<EdgeKey2D, 3> local_edges{{
                canonical_edge(element[0], element[1]),
                canonical_edge(element[1], element[2]),
                canonical_edge(element[2], element[0])
            }};

            for (std::size_t local_edge = 0; local_edge < 3; ++local_edge) {
                const auto iterator = edge_indices.find(local_edges[local_edge]);
                if (iterator == edge_indices.end()) {
                    throw std::logic_error("Internal H1 edge numbering inconsistency.");
                }
                mapping.dofs[3U + local_edge] =
                    vertex_count_ + iterator->second;
            }
        }

        element_dofs_.push_back(mapping);
    }
}

std::size_t H1DofManager2D::vertex_dof(
    const std::size_t vertex
) const {
    if (vertex >= vertex_count_) {
        throw std::out_of_range("H1 vertex DOF references a missing mesh node.");
    }
    return vertex;
}

std::optional<std::size_t> H1DofManager2D::edge_dof(
    const std::size_t first,
    const std::size_t second
) const {
    if (order_ == H1TriangleOrder::P1) {
        return std::nullopt;
    }

    if (first >= vertex_count_ || second >= vertex_count_) {
        throw std::out_of_range("H1 edge references a missing mesh node.");
    }
    if (first == second) {
        throw std::invalid_argument("An H1 edge must have two distinct nodes.");
    }

    const EdgeKey2D key = canonical_edge(first, second);

    // edge_keys_ is stored in discovery order, so use an explicit linear
    // lookup. The public API is independent of the backing index structure.
    for (std::size_t index = 0; index < edge_keys_.size(); ++index) {
        if (edge_keys_[index] == key) {
            return vertex_count_ + index;
        }
    }

    return std::nullopt;
}

H1DofEntity H1DofManager2D::dof_entity(const std::size_t dof) const {
    if (dof >= dof_count_) {
        throw std::out_of_range("H1 DOF index is out of range.");
    }
    return dof < vertex_count_ ? H1DofEntity::Vertex : H1DofEntity::Edge;
}

const H1ElementDofMap2D& H1DofManager2D::element_dofs(
    const std::size_t element
) const {
    if (element >= element_dofs_.size()) {
        throw std::out_of_range("H1 element index is out of range.");
    }
    return element_dofs_[element];
}

H1FunctionSpace2D::H1FunctionSpace2D(
    const TriangleMesh2D& mesh,
    const H1TriangleOrder order
)
    : mesh_(&mesh),
      dof_manager_(mesh, order) {}

std::vector<double> H1FunctionSpace2D::interpolate(
    const std::function<double(double, double)>& field
) const {
    if (!field) {
        throw std::invalid_argument("An H1 interpolation field must be provided.");
    }

    std::vector<double> values(dimension(), 0.0);

    for (std::size_t vertex = 0; vertex < mesh_->nodes.size(); ++vertex) {
        const Point2D& point = mesh_->nodes[vertex];
        values[vertex] = field(point.x, point.y);
    }

    if (order() == H1TriangleOrder::P2) {
        for (std::size_t edge = 0; edge < dof_manager_.edge_keys().size(); ++edge) {
            const EdgeKey2D& key = dof_manager_.edge_keys()[edge];
            const Point2D& p0 = mesh_->nodes[key.first];
            const Point2D& p1 = mesh_->nodes[key.second];
            const double x = 0.5 * (p0.x + p1.x);
            const double y = 0.5 * (p0.y + p1.y);
            values[mesh_->nodes.size() + edge] = field(x, y);
        }
    }

    return values;
}

} // namespace nova::fem
