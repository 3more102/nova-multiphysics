#include "nova/fem/error_norms_2d.hpp"

#include "nova/fem/triangle_p1.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace nova::fem {
namespace {

struct GaussPoint1D {
    double coordinate;
    double weight;
};

// Five-point Gauss-Legendre rule mapped to [0,1], exact through degree 9.
constexpr std::array<GaussPoint1D, 5> gauss_legendre_5{{
    {0.046910077030668, 0.118463442528095},
    {0.230765344947158, 0.239314335249683},
    {0.5,               0.284444444444444},
    {0.769234655052842, 0.239314335249683},
    {0.953089922969332, 0.118463442528095}
}};

} // namespace

ErrorNorms2D compute_error_norms_2d(
    const TriangleMesh2D& mesh,
    const std::span<const double> nodal_values,
    const ExactScalarField2D& exact_solution
) {
    mesh.validate();

    if (nodal_values.size() != mesh.nodes.size()) {
        throw std::invalid_argument("Error-norm field size must match the mesh node count.");
    }
    if (!exact_solution.value ||
        !exact_solution.derivative_x ||
        !exact_solution.derivative_y) {
        throw std::invalid_argument("Exact solution and both spatial derivatives are required.");
    }

    for (const double value : nodal_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Numerical field values must be finite.");
        }
    }

    double l2_squared = 0.0;
    double h1_seminorm_squared = 0.0;

    for (const Triangle2D& element : mesh.elements) {
        const Point2D& p0 = mesh.nodes[element[0]];
        const Point2D& p1 = mesh.nodes[element[1]];
        const Point2D& p2 = mesh.nodes[element[2]];

        const double determinant =
            (p1.x - p0.x) * (p2.y - p0.y) -
            (p2.x - p0.x) * (p1.y - p0.y);

        const std::array<std::array<double, 2>, 3> gradients{{
            {{(p1.y - p2.y) / determinant, (p2.x - p1.x) / determinant}},
            {{(p2.y - p0.y) / determinant, (p0.x - p2.x) / determinant}},
            {{(p0.y - p1.y) / determinant, (p1.x - p0.x) / determinant}}
        }};

        const std::array<double, 3> local_values{{
            nodal_values[element[0]],
            nodal_values[element[1]],
            nodal_values[element[2]]
        }};

        double numerical_gradient_x = 0.0;
        double numerical_gradient_y = 0.0;
        for (std::size_t local = 0; local < 3; ++local) {
            numerical_gradient_x += local_values[local] * gradients[local][0];
            numerical_gradient_y += local_values[local] * gradients[local][1];
        }

        for (const GaussPoint1D& point_u : gauss_legendre_5) {
            for (const GaussPoint1D& point_v : gauss_legendre_5) {
                // Duffy map from the unit square to the reference triangle.
                const double xi = point_u.coordinate;
                const double eta = (1.0 - xi) * point_v.coordinate;
                const double reference_weight =
                    point_u.weight * point_v.weight * (1.0 - xi);

                const std::array<double, 3> shape =
                    triangle_p1_shape_values(xi, eta);

                const double x =
                    shape[0] * p0.x + shape[1] * p1.x + shape[2] * p2.x;
                const double y =
                    shape[0] * p0.y + shape[1] * p1.y + shape[2] * p2.y;

                const double numerical_value =
                    shape[0] * local_values[0] +
                    shape[1] * local_values[1] +
                    shape[2] * local_values[2];

                const double exact_value = exact_solution.value(x, y);
                const double exact_gradient_x = exact_solution.derivative_x(x, y);
                const double exact_gradient_y = exact_solution.derivative_y(x, y);

                if (!std::isfinite(exact_value) ||
                    !std::isfinite(exact_gradient_x) ||
                    !std::isfinite(exact_gradient_y)) {
                    throw std::invalid_argument("Exact solution returned a non-finite value.");
                }

                const double value_error = numerical_value - exact_value;
                const double gradient_error_x =
                    numerical_gradient_x - exact_gradient_x;
                const double gradient_error_y =
                    numerical_gradient_y - exact_gradient_y;
                const double physical_weight =
                    std::abs(determinant) * reference_weight;

                l2_squared += value_error * value_error * physical_weight;
                h1_seminorm_squared +=
                    (gradient_error_x * gradient_error_x +
                     gradient_error_y * gradient_error_y) * physical_weight;

                if (!std::isfinite(l2_squared) ||
                    !std::isfinite(h1_seminorm_squared)) {
                    throw std::overflow_error("Error-norm integration overflowed.");
                }
            }
        }
    }

    const double h1_squared = l2_squared + h1_seminorm_squared;
    if (!std::isfinite(h1_squared)) {
        throw std::overflow_error("Combined H1 error-norm integration overflowed.");
    }

    const double l2 = std::sqrt(l2_squared);
    const double h1_seminorm = std::sqrt(h1_seminorm_squared);
    return {l2, h1_seminorm, std::sqrt(h1_squared)};
}

} // namespace nova::fem
