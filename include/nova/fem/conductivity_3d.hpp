#pragma once

#include <functional>

namespace nova::fem {

using SpatialScalarField3D = std::function<double(double, double, double)>;

// Components in the global Cartesian frame. Each off-diagonal entry represents
// both symmetric positions. Conductivity must be positive definite at each
// integration point. The public problem APIs accept either a scalar or a tensor.
struct SymmetricTensor3D {
    double xx{}, yy{}, zz{};
    double xy{}, xz{}, yz{};
};

using SymmetricTensorField3D = std::function<SymmetricTensor3D(double, double, double)>;

} // namespace nova::fem
