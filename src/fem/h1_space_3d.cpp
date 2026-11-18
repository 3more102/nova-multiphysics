#include "nova/fem/h1_space_3d.hpp"

#include <cmath>
#include <stdexcept>

namespace nova::fem {

H1DofManager3D::H1DofManager3D(const TetrahedralMesh3D& mesh)
    : dof_count_(mesh.nodes.size()) {
    mesh.validate();
    element_dofs_.reserve(mesh.elements.size());
    for (const Tetrahedron3D& tetrahedron : mesh.elements) {
        element_dofs_.push_back(tetrahedron);
    }
}

const Tetrahedron3D& H1DofManager3D::element_dofs(
    const std::size_t element
) const {
    return element_dofs_.at(element);
}

H1FunctionSpace3D::H1FunctionSpace3D(const TetrahedralMesh3D& mesh)
    : mesh_(&mesh), dof_manager_(mesh) {}

std::vector<double> H1FunctionSpace3D::interpolate(
    const std::function<double(double, double, double)>& field
) const {
    if (!field) {
        throw std::invalid_argument("H1 3D interpolation requires a field callback.");
    }
    std::vector<double> values;
    values.reserve(mesh_->nodes.size());
    for (const Point3D& node : mesh_->nodes) {
        const double value = field(node.x, node.y, node.z);
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "H1 3D interpolation field returned a non-finite value."
            );
        }
        values.push_back(value);
    }
    return values;
}

} // namespace nova::fem
