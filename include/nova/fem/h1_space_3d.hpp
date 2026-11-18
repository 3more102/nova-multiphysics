#pragma once

#include "nova/fem/mesh3d.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

// Continuous H1 P1 numbering for tetrahedral meshes. The initial 3D space
// uses one global degree of freedom per mesh vertex.
class H1DofManager3D {
public:
    explicit H1DofManager3D(const TetrahedralMesh3D& mesh);

    [[nodiscard]] std::size_t dof_count() const noexcept { return dof_count_; }
    [[nodiscard]] std::size_t element_count() const noexcept {
        return element_dofs_.size();
    }
    [[nodiscard]] const Tetrahedron3D& element_dofs(
        std::size_t element
    ) const;

private:
    std::size_t dof_count_{};
    std::vector<Tetrahedron3D> element_dofs_;
};

class H1FunctionSpace3D {
public:
    explicit H1FunctionSpace3D(const TetrahedralMesh3D& mesh);

    [[nodiscard]] const TetrahedralMesh3D& mesh() const noexcept {
        return *mesh_;
    }
    [[nodiscard]] std::size_t dimension() const noexcept {
        return dof_manager_.dof_count();
    }
    [[nodiscard]] const H1DofManager3D& dofs() const noexcept {
        return dof_manager_;
    }

    [[nodiscard]] std::vector<double> interpolate(
        const std::function<double(double, double, double)>& field
    ) const;

private:
    const TetrahedralMesh3D* mesh_{};
    H1DofManager3D dof_manager_;
};

} // namespace nova::fem
