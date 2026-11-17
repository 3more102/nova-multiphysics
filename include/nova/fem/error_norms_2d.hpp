#pragma once

#include "nova/fem/mesh2d.hpp"

#include <functional>
#include <span>

namespace nova::fem {

struct ExactScalarField2D {
    std::function<double(double, double)> value;
    std::function<double(double, double)> derivative_x;
    std::function<double(double, double)> derivative_y;
};

struct ErrorNorms2D {
    double l2{0.0};
    double h1_seminorm{0.0};
    double h1{0.0};
};

// Computes L2, H1-seminorm, and H1 errors for a nodal P1 field.
// Uses a 5x5 Gauss-Legendre rule with a Duffy map to integrate each triangle.
[[nodiscard]] ErrorNorms2D compute_error_norms_2d(
    const TriangleMesh2D& mesh,
    std::span<const double> nodal_values,
    const ExactScalarField2D& exact_solution
);

} // namespace nova::fem
