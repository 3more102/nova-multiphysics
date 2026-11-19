#pragma once

#include "nova/fem/mesh3d.hpp"

#include <functional>
#include <span>

namespace nova::fem {

struct ExactScalarField3D {
    std::function<double(double, double, double)> value;
    std::function<double(double, double, double)> derivative_x;
    std::function<double(double, double, double)> derivative_y;
    std::function<double(double, double, double)> derivative_z;
};

struct ErrorNorms3D {
    double l2{0.0};
    double h1_seminorm{0.0};
    double h1{0.0};
};

// Verification-only L2, H1-seminorm, and H1 errors for a nodal P1
// tetrahedral field. Uses 5^3 Gauss-Legendre samples through a 3D
// Duffy map, accurate for polynomial error integrands through degree 4.
// Requires an exact value and exact Cartesian gradient.
[[nodiscard]] ErrorNorms3D compute_error_norms_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const double> nodal_values,
    const ExactScalarField3D& exact_solution
);

} // namespace nova::fem
