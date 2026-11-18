#include "nova/fem/residual_error_estimator_2d.hpp"

#include "nova/fem/triangle_p1.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using EdgeKey = std::pair<std::size_t, std::size_t>;

struct ElementGeometry2D {
    double j00{};
    double j01{};
    double j10{};
    double j11{};
    double determinant{};
    double absolute_determinant{};
    double diameter{};
};

struct Gradient2D {
    double x{};
    double y{};
};

struct EdgeUse2D {
    std::size_t element{};
    std::size_t first{};
    std::size_t second{};
};

struct EdgeQuadraturePoint {
    double t{};
    double weight{};
};

[[nodiscard]] EdgeKey canonical_edge(
    const std::size_t first, const std::size_t second
) {
    return std::minmax(first, second);
}

[[nodiscard]] double distance(
    const Point2D& first, const Point2D& second
) {
    return std::hypot(second.x - first.x, second.y - first.y);
}

[[nodiscard]] ElementGeometry2D element_geometry(
    const TriangleMesh2D& mesh, const std::size_t element
) {
    const Triangle2D& triangle = mesh.elements[element];
    const Point2D& p0 = mesh.nodes[triangle[0]];
    const Point2D& p1 = mesh.nodes[triangle[1]];
    const Point2D& p2 = mesh.nodes[triangle[2]];

    ElementGeometry2D geometry;
    geometry.j00 = p1.x - p0.x;
    geometry.j01 = p2.x - p0.x;
    geometry.j10 = p1.y - p0.y;
    geometry.j11 = p2.y - p0.y;
    geometry.determinant =
        geometry.j00 * geometry.j11 - geometry.j01 * geometry.j10;
    geometry.absolute_determinant = std::abs(geometry.determinant);
    geometry.diameter = std::max({
        distance(p0, p1), distance(p1, p2), distance(p2, p0)
    });
    if (!std::isfinite(geometry.determinant) ||
        geometry.determinant == 0.0 ||
        !std::isfinite(geometry.diameter) ||
        geometry.diameter <= 0.0) {
        throw std::invalid_argument(
            "Residual estimator encountered degenerate element geometry."
        );
    }
    return geometry;
}

[[nodiscard]] std::array<double, 2> reference_coordinates(
    const TriangleMesh2D& mesh,
    const ElementGeometry2D& geometry,
    const std::size_t element,
    const Point2D& point
) {
    const Triangle2D& triangle = mesh.elements[element];
    const Point2D& p0 = mesh.nodes[triangle[0]];
    const double dx = point.x - p0.x;
    const double dy = point.y - p0.y;
    const double r = (
        dx * geometry.j11 - dy * geometry.j01
    ) / geometry.determinant;
    const double s = (
        geometry.j00 * dy - geometry.j10 * dx
    ) / geometry.determinant;
    constexpr double tolerance = 1.0e-10;
    if (!std::isfinite(r) || !std::isfinite(s) ||
        r < -tolerance || s < -tolerance ||
        r + s > 1.0 + tolerance) {
        throw std::invalid_argument(
            "Interior edge quadrature point is outside an adjacent element."
        );
    }
    return {r, s};
}

[[nodiscard]] Gradient2D evaluate_gradient(
    const TriangleMesh2D& mesh,
    const H1DofManager2D& dofs,
    const H1TriangleOrder order,
    const std::span<const double> values,
    const std::size_t element,
    const Point2D& point
) {
    const ElementGeometry2D geometry = element_geometry(mesh, element);
    const auto reference = reference_coordinates(
        mesh, geometry, element, point
    );
    const TriangleH1BasisEvaluation basis =
        evaluate_triangle_h1_basis(order, reference[0], reference[1]);
    const H1ElementDofMap2D& local = dofs.element_dofs(element);

    double gradient_r = 0.0;
    double gradient_s = 0.0;
    for (std::size_t i = 0; i < basis.count; ++i) {
        gradient_r += values[local[i]] * basis.reference_gradients[i][0];
        gradient_s += values[local[i]] * basis.reference_gradients[i][1];
    }

    const double gradient_x = (
        geometry.j11 * gradient_r - geometry.j10 * gradient_s
    ) / geometry.determinant;
    const double gradient_y = (
        -geometry.j01 * gradient_r + geometry.j00 * gradient_s
    ) / geometry.determinant;
    if (!std::isfinite(gradient_x) || !std::isfinite(gradient_y)) {
        throw std::overflow_error("Residual estimator gradient is non-finite.");
    }
    return {gradient_x, gradient_y};
}

[[nodiscard]] double evaluate_element_laplacian(
    const TriangleMesh2D& mesh,
    const H1DofManager2D& dofs,
    const H1TriangleOrder order,
    const std::span<const double> values,
    const std::size_t element
) {
    if (order == H1TriangleOrder::P1) {
        return 0.0;
    }

    // Reference Hessians [drr, drs, dss] for the six local P2 basis
    // functions in the same vertex/edge ordering as triangle_h1.hpp.
    constexpr std::array<double, 6> hessian_rr{{
        4.0, 4.0, 0.0, -8.0, 0.0, 0.0
    }};
    constexpr std::array<double, 6> hessian_rs{{
        4.0, 0.0, 0.0, -4.0, 4.0, -4.0
    }};
    constexpr std::array<double, 6> hessian_ss{{
        4.0, 0.0, 4.0, 0.0, 0.0, -8.0
    }};

    const ElementGeometry2D geometry = element_geometry(mesh, element);
    const H1ElementDofMap2D& local = dofs.element_dofs(element);
    double h_rr = 0.0;
    double h_rs = 0.0;
    double h_ss = 0.0;
    for (std::size_t i = 0; i < 6U; ++i) {
        const double value = values[local[i]];
        h_rr += value * hessian_rr[i];
        h_rs += value * hessian_rs[i];
        h_ss += value * hessian_ss[i];
    }

    // For an affine element map, H_x = J^{-T} H_ref J^{-1}.
    const double r_x = geometry.j11 / geometry.determinant;
    const double r_y = -geometry.j01 / geometry.determinant;
    const double s_x = -geometry.j10 / geometry.determinant;
    const double s_y = geometry.j00 / geometry.determinant;
    const double laplacian =
        (r_x * r_x + r_y * r_y) * h_rr +
        2.0 * (r_x * s_x + r_y * s_y) * h_rs +
        (s_x * s_x + s_y * s_y) * h_ss;
    if (!std::isfinite(laplacian)) {
        throw std::overflow_error("P2 element Laplacian is non-finite.");
    }
    return laplacian;
}

[[nodiscard]] std::array<double, 2> outward_normal(
    const TriangleMesh2D& mesh,
    const EdgeUse2D& use
) {
    const ElementGeometry2D geometry = element_geometry(mesh, use.element);
    const Point2D& first = mesh.nodes[use.first];
    const Point2D& second = mesh.nodes[use.second];
    const double dx = second.x - first.x;
    const double dy = second.y - first.y;
    const double length = std::hypot(dx, dy);
    if (!std::isfinite(length) || length <= 0.0) {
        throw std::invalid_argument("Interior edge has invalid length.");
    }
    const double orientation = geometry.determinant > 0.0 ? 1.0 : -1.0;
    return {
        orientation * dy / length,
        -orientation * dx / length
    };
}

[[nodiscard]] constexpr std::array<EdgeQuadraturePoint, 3>
edge_quadrature() {
    const double offset = 0.5 * std::sqrt(3.0 / 5.0);
    return {{
        {0.5 - offset, 5.0 / 18.0},
        {0.5, 8.0 / 18.0},
        {0.5 + offset, 5.0 / 18.0}
    }};
}

} // namespace

DiffusionResidualIndicators2D estimate_diffusion_residual_error_2d(
    const TriangleMesh2D& mesh,
    const H1TriangleOrder order,
    const std::span<const double> dof_values,
    const std::span<const double> element_conductivity,
    const std::function<double(double, double)>& source
) {
    mesh.validate();
    const H1FunctionSpace2D space(mesh, order);
    const H1DofManager2D& dofs = space.dofs();
    if (dof_values.size() != space.dimension()) {
        throw std::invalid_argument(
            "Residual estimator received the wrong H1 DOF count."
        );
    }
    if (element_conductivity.size() != mesh.elements.size()) {
        throw std::invalid_argument(
            "Residual estimator conductivity count must match the element count."
        );
    }
    if (!source) {
        throw std::invalid_argument("Residual estimator requires a source function.");
    }
    for (const double value : dof_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "Residual estimator state must contain finite values."
            );
        }
    }
    for (const double conductivity : element_conductivity) {
        if (!std::isfinite(conductivity) || conductivity <= 0.0) {
            throw std::invalid_argument(
                "Element conductivity must be finite and strictly positive."
            );
        }
    }

    std::vector<double> indicator_squared(mesh.elements.size(), 0.0);
    std::vector<ElementGeometry2D> geometries;
    geometries.reserve(mesh.elements.size());

    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const ElementGeometry2D geometry = element_geometry(mesh, element);
        geometries.push_back(geometry);
        const double conductivity = element_conductivity[element];
        const double laplacian = evaluate_element_laplacian(
            mesh, dofs, order, dof_values, element
        );

        double residual_integral = 0.0;
        const Triangle2D& triangle = mesh.elements[element];
        const Point2D& p0 = mesh.nodes[triangle[0]];
        for (const TriangleQuadraturePoint& quadrature :
             triangle_degree_four_quadrature) {
            const double x = p0.x +
                geometry.j00 * quadrature.r + geometry.j01 * quadrature.s;
            const double y = p0.y +
                geometry.j10 * quadrature.r + geometry.j11 * quadrature.s;
            const double source_value = source(x, y);
            if (!std::isfinite(source_value)) {
                throw std::invalid_argument(
                    "Residual estimator source must be finite."
                );
            }
            const double residual = source_value + conductivity * laplacian;
            const double contribution =
                residual * residual * geometry.absolute_determinant *
                quadrature.weight;
            residual_integral += contribution;
            if (!std::isfinite(residual_integral)) {
                throw std::overflow_error(
                    "Residual estimator cell integral overflowed."
                );
            }
        }

        indicator_squared[element] =
            geometry.diameter * geometry.diameter * residual_integral;
        if (!std::isfinite(indicator_squared[element])) {
            throw std::overflow_error(
                "Residual estimator cell contribution overflowed."
            );
        }
    }

    std::map<EdgeKey, std::vector<EdgeUse2D>> edge_uses;
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const Triangle2D& triangle = mesh.elements[element];
        edge_uses[canonical_edge(triangle[0], triangle[1])].push_back(
            {element, triangle[0], triangle[1]}
        );
        edge_uses[canonical_edge(triangle[1], triangle[2])].push_back(
            {element, triangle[1], triangle[2]}
        );
        edge_uses[canonical_edge(triangle[2], triangle[0])].push_back(
            {element, triangle[2], triangle[0]}
        );
    }

    constexpr auto quadrature = edge_quadrature();
    for (const auto& [edge, uses] : edge_uses) {
        if (uses.size() == 1U) {
            // Natural-boundary residual terms are intentionally omitted.
            continue;
        }
        if (uses.size() != 2U) {
            throw std::invalid_argument(
                "Residual estimator requires a manifold triangular mesh."
            );
        }

        const EdgeUse2D& first_use = uses[0];
        const EdgeUse2D& second_use = uses[1];
        const Point2D& first_point = mesh.nodes[edge.first];
        const Point2D& second_point = mesh.nodes[edge.second];
        const double edge_length = distance(first_point, second_point);
        if (!std::isfinite(edge_length) || edge_length <= 0.0) {
            throw std::invalid_argument(
                "Residual estimator encountered an invalid edge length."
            );
        }
        const auto first_normal = outward_normal(mesh, first_use);
        const auto second_normal = outward_normal(mesh, second_use);
        double jump_integral = 0.0;

        for (const EdgeQuadraturePoint& quadrature_point : quadrature) {
            const Point2D point{
                (1.0 - quadrature_point.t) * first_point.x +
                    quadrature_point.t * second_point.x,
                (1.0 - quadrature_point.t) * first_point.y +
                    quadrature_point.t * second_point.y
            };
            const Gradient2D first_gradient = evaluate_gradient(
                mesh, dofs, order, dof_values, first_use.element, point
            );
            const Gradient2D second_gradient = evaluate_gradient(
                mesh, dofs, order, dof_values, second_use.element, point
            );
            const double first_flux = element_conductivity[first_use.element] *
                (first_gradient.x * first_normal[0] +
                 first_gradient.y * first_normal[1]);
            const double second_flux = element_conductivity[second_use.element] *
                (second_gradient.x * second_normal[0] +
                 second_gradient.y * second_normal[1]);
            const double jump = first_flux + second_flux;
            jump_integral += jump * jump * edge_length * quadrature_point.weight;
            if (!std::isfinite(jump_integral)) {
                throw std::overflow_error(
                    "Residual estimator edge-jump integral overflowed."
                );
            }
        }

        const double edge_contribution = 0.5 * edge_length * jump_integral;
        indicator_squared[first_use.element] += edge_contribution;
        indicator_squared[second_use.element] += edge_contribution;
        if (!std::isfinite(indicator_squared[first_use.element]) ||
            !std::isfinite(indicator_squared[second_use.element])) {
            throw std::overflow_error(
                "Residual estimator edge contribution overflowed."
            );
        }
    }

    DiffusionResidualIndicators2D result;
    result.element_errors.resize(mesh.elements.size(), 0.0);
    double total_squared = 0.0;
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        result.element_errors[element] = std::sqrt(indicator_squared[element]);
        total_squared += indicator_squared[element];
        if (!std::isfinite(total_squared)) {
            throw std::overflow_error(
                "Residual estimator global norm overflowed."
            );
        }
    }
    result.total_error = std::sqrt(total_squared);
    return result;
}

} // namespace nova::fem
