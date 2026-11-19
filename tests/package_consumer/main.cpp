#include "nova/fem/poisson_1d.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <cmath>
#include <vector>

int main() {
    // Exercises FEM and sparse algebra across an installed library boundary.
    nova::fem::Poisson1DProblem p;
    p.elements = 16;
    p.source = [](double) { return 2.0; };
    const auto solution = nova::fem::solve_poisson_1d(p);
    if (solution.values.size() != 17 ||
        std::abs(solution.values[8] - 0.25) > 1.0e-7) {
        return 1;
    }

    const auto a = nova::linalg::CsrMatrix::from_triplets(2, 2,
        {{0, 0, 2.0}, {1, 1, 3.0}});
    const std::vector<double> x{4.0, 5.0};
    const auto y = a.multiply(x);
    return y.size() == 2 && std::abs(y[0] - 8.0) < 1.0e-12 &&
        std::abs(y[1] - 15.0) < 1.0e-12 ? 0 : 2;
}
