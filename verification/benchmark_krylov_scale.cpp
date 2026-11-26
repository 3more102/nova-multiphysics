// Standalone driver: compile this same file against each library under comparison.
#include "nova/linalg/bicgstab.hpp"
#include "nova/linalg/gmres.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <span>
#include <string_view>
#include <vector>

namespace {
constexpr std::size_t dimension = 1024;
constexpr double diagonal = 3.5, lower = -1.25, upper = -0.5;
constexpr double rtol = 1e-11;
constexpr int batch_size = 64;
constexpr int samples = 7;

double stencil(std::span<const double> x, std::size_t i) {
    const double left = i == 0 ? 0.0 : x[i - 1];
    const double right = i + 1 == x.size() ? 0.0 : x[i + 1];
    return diagonal * x[i] + lower * left + upper * right;
}

template<class Solve>
bool measure(std::string_view name, const std::vector<double>& exact,
             const std::vector<double>& unit_load, Solve solve) {
    bool verified = true;
    for (const int exponent : {0, -660, 660}) {
        const double amplitude = std::ldexp(1.0, exponent);
        auto rhs = unit_load;
        for (double& v : rhs) v *= amplitude;

        std::cout << name << " source_exp=" << exponent;
        try {
            auto result = solve(rhs); // warmup for timings, single solve for the scale checks
            if (exponent == 0 && result.converged) {
                std::array<double, samples> times{};
                for (double& sample : times) {
                    const auto start = std::chrono::steady_clock::now();
                    for (int rep = 0; rep < batch_size; ++rep) result = solve(rhs);
                    const auto elapsed = std::chrono::steady_clock::now() - start;
                    sample = std::chrono::duration<double, std::micro>(elapsed).count() / batch_size;
                }
                std::sort(times.begin(), times.end());
                std::cout << " median_us=" << times[times.size() / 2];
            }
            // Normalize first, then apply the independent stencil in physical row order.
            // Every entry of the unit solution and unit load is exactly representable.
            auto x = result.solution;
            for (double& v : x) v /= amplitude;
            double error = 0, solution_norm = 0, residual = 0, load_norm = 0;
            for (std::size_t i = 0; i < dimension; ++i) {
                error = std::hypot(error, x.at(i) - exact[i]);
                solution_norm = std::hypot(solution_norm, exact[i]);
                residual = std::hypot(residual, unit_load[i] - stencil(x, i));
                load_norm = std::hypot(load_norm, unit_load[i]);
            }
            const double rel_error = error / solution_norm;
            const double rel_residual = residual / load_norm;
            const bool accurate = result.converged && rel_error <= 4 * rtol && rel_residual <= rtol;
            std::cout << " iterations=" << result.iterations << " solver_converged=" << result.converged
                      << " rel_error=" << rel_error << " rel_residual=" << rel_residual
                      << " verified=" << accurate;
            if (!result.converged) std::cout << " diagnostic=" << result.message;
            verified &= accurate;
        } catch (const std::exception& e) {
            std::cout << " verified=0 exception=" << e.what();
            verified = false;
        }
        std::cout << '\n';
    }
    return verified;
}
} // namespace

int main() {
    using namespace nova::linalg;
    std::vector<double> exact(dimension), unit_load(dimension);
    std::vector<Triplet> entries;
    for (std::size_t i = 0; i < dimension; ++i) {
        exact[i] = (static_cast<int>((37 * i) % 31) - 15) / 16.0;
        if (i > 0) entries.push_back({i, i - 1, lower});
        entries.push_back({i, i, diagonal});
        if (i + 1 < dimension) entries.push_back({i, i + 1, upper});
    }
    for (std::size_t i = 0; i < dimension; ++i) unit_load[i] = stencil(exact, i);
    const auto matrix = CsrMatrix::from_triplets(dimension, dimension, entries);

    GmresOptions options;
    options.restart = 16;
    options.solver.absolute_tolerance = 0.0;
    options.solver.relative_tolerance = rtol;
    std::cout << std::setprecision(7) << "dofs=" << dimension << " nnz=" << matrix.non_zeros()
              << " restart=" << options.restart << " rtol=" << rtol << " atol=0 preconditioner=none"
              << " samples=" << samples << " solves_per_sample=" << batch_size << '\n';
    const bool gmres = measure("GMRES", exact, unit_load, [&](std::span<const double> rhs) {
        return solve_gmres(matrix, rhs, options);
    });
    const bool bicgstab = measure("BiCGSTAB", exact, unit_load, [&](std::span<const double> rhs) {
        return solve_bicgstab(matrix, rhs, options.solver);
    });
    return gmres && bicgstab ? 0 : 1;
}
