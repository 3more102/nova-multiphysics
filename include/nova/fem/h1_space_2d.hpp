#pragma once

#include "nova/fem/mesh2d.hpp"
#include "nova/fem/triangle_h1.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

namespace nova::fem {

using EdgeKey2D = std::pair<std::size_t, std::size_t>;

enum class H1DofEntity {
    Vertex,
    Edge
};

struct H1ElementDofMap2D {
    std::array<std::size_t, 6> dofs{};
    std::size_t count{};

    [[nodiscard]] std::size_t operator[](const std::size_t local_index) const {
        return dofs.at(local_index);
    }
};

// Global continuous H1 Lagrange DOF numbering for a triangular 2D mesh.
//
// P1:
//   global DOFs = mesh vertices.
//
// P2:
//   global DOFs = mesh vertices followed by unique mesh edges.
//   Edge DOFs are shared by adjacent elements through canonicalized
//   (min(node_i,node_j), max(node_i,node_j)) edge keys.
class H1DofManager2D {
public:
    H1DofManager2D(
        const TriangleMesh2D& mesh,
        H1TriangleOrder order
    );

    [[nodiscard]] H1TriangleOrder order() const noexcept { return order_; }
    [[nodiscard]] std::size_t dof_count() const noexcept { return dof_count_; }
    [[nodiscard]] std::size_t element_count() const noexcept {
        return element_dofs_.size();
    }
    [[nodiscard]] std::size_t element_dof_count() const noexcept {
        return triangle_h1_dof_count(order_);
    }
    [[nodiscard]] std::size_t vertex_dof_count() const noexcept {
        return vertex_count_;
    }
    [[nodiscard]] std::size_t edge_dof_count() const noexcept {
        return order_ == H1TriangleOrder::P2 ? edge_keys_.size() : 0U;
    }

    [[nodiscard]] std::size_t vertex_dof(std::size_t vertex) const;
    [[nodiscard]] std::optional<std::size_t> edge_dof(
        std::size_t first,
        std::size_t second
    ) const;

    [[nodiscard]] H1DofEntity dof_entity(std::size_t dof) const;

    [[nodiscard]] const H1ElementDofMap2D& element_dofs(
        std::size_t element
    ) const;

    [[nodiscard]] const std::vector<EdgeKey2D>& edge_keys() const noexcept {
        return edge_keys_;
    }

private:
    H1TriangleOrder order_{H1TriangleOrder::P1};
    std::size_t vertex_count_{};
    std::size_t dof_count_{};
    std::vector<EdgeKey2D> edge_keys_;
    std::vector<H1ElementDofMap2D> element_dofs_;
};

class H1FunctionSpace2D {
public:
    H1FunctionSpace2D(
        const TriangleMesh2D& mesh,
        H1TriangleOrder order
    );

    [[nodiscard]] const TriangleMesh2D& mesh() const noexcept { return *mesh_; }
    [[nodiscard]] H1TriangleOrder order() const noexcept {
        return dof_manager_.order();
    }
    [[nodiscard]] std::size_t dimension() const noexcept {
        return dof_manager_.dof_count();
    }
    [[nodiscard]] std::size_t element_dof_count() const noexcept {
        return dof_manager_.element_dof_count();
    }
    [[nodiscard]] const H1DofManager2D& dofs() const noexcept {
        return dof_manager_;
    }

    // Nodal interpolation. P1 samples vertices; P2 samples vertices and
    // edge midpoints, matching the global nodal DOF layout.
    [[nodiscard]] std::vector<double> interpolate(
        const std::function<double(double, double)>& field
    ) const;

private:
    const TriangleMesh2D* mesh_{};
    H1DofManager2D dof_manager_;
};

} // namespace nova::fem
