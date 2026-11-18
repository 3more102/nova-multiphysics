#include "nova/fem/error_indicators_2d.hpp"

#include "nova/fem/h1_space_2d.hpp"
#include "nova/fem/triangle_h1.hpp"
#include "nova/fem/triangle_p1.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

struct Gradient2D {
    double x{};
    double y{};
};

[[nodiscard]] Gradient2D evaluate_element_gradient(
    const TriangleMesh2D& mesh,
    const H1DofManager2D& dofs,
    const H1TriangleOrder order,
    const std::span<const double> values,
    const std::size_t element,
    const double r,
    const double s
) {
    const Triangle2D& triangle = mesh.elements[element];
    const Point2D& p0 = mesh.nodes[triangle[0]];
    const Point2D& p1 = mesh.nodes[triangle[1]];
    const Point2D& p2 = mesh.nodes[triangle[2]];
    const double j00 = p1.x - p0.x;
    const double j01 = p2.x - p0.x;
    const double j10 = p1.y - p0.y;
    const double j11 = p2.y - p0.y;
    const double determinant = j00 * j11 - j01 * j10;
    if (!std::isfinite(determinant) || determinant == 0.0) {
        throw std::invalid_argument("Gradient recovery encountered a degenerate element.");
    }

    const TriangleH1BasisEvaluation basis =
        evaluate_triangle_h1_basis(order, r, s);
    const H1ElementDofMap2D& local_dofs = dofs.element_dofs(element);
    Gradient2D gradient;
    for (std::size_t local = 0; local < basis.count; ++local) {
        const double dr = basis.reference_gradients[local][0];
        const double ds = basis.reference_gradients[local][1];
        const double dx = (j11 * dr - j10 * ds) / determinant;
        const double dy = (-j01 * dr + j00 * ds) / determinant;
        gradient.x += values[local_dofs[local]] * dx;
        gradient.y += values[local_dofs[local]] * dy;
    }
    if (!std::isfinite(gradient.x) || !std::isfinite(gradient.y)) {
        throw std::overflow_error("Finite-element gradient is non-finite.");
    }
    return gradient;
}

[[nodiscard]] double element_area(
    const TriangleMesh2D& mesh,
    const std::size_t element
) {
    const Triangle2D& triangle = mesh.elements[element];
    const Point2D& a = mesh.nodes[triangle[0]];
    const Point2D& b = mesh.nodes[triangle[1]];
    const Point2D& c = mesh.nodes[triangle[2]];
    return 0.5 * std::abs(
        (b.x - a.x) * (c.y - a.y) -
        (b.y - a.y) * (c.x - a.x)
    );
}

} // namespace

H1GradientRecoveryIndicators2D
estimate_h1_gradient_recovery_error_2d(
    const TriangleMesh2D& mesh,
    const H1TriangleOrder order,
    const std::span<const double> dof_values,
    std::function<double(double, double)> conductivity
) {
    mesh.validate();
    const H1FunctionSpace2D space(mesh, order);
    const H1DofManager2D& dofs = space.dofs();
    if (dof_values.size() != space.dimension()) {
        throw std::invalid_argument("Gradient recovery state has the wrong H1 DOF count.");
    }
    for (const double value : dof_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Gradient recovery state must be finite.");
        }
    }

    std::vector<Gradient2D> accumulated(mesh.nodes.size());
    std::vector<double> accumulated_weights(mesh.nodes.size(), 0.0);
    constexpr std::array<std::array<double, 2>, 3> vertex_reference{{
        {{0.0, 0.0}}, {{1.0, 0.0}}, {{0.0, 1.0}}
    }};

    // Recover a continuous nodal gradient by area-weighted patch averaging.
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const double area = element_area(mesh, element);
        const Triangle2D& triangle = mesh.elements[element];
        for (std::size_t local_vertex = 0; local_vertex < 3U; ++local_vertex) {
            const Gradient2D gradient = evaluate_element_gradient(
                mesh, dofs, order, dof_values, element,
                vertex_reference[local_vertex][0],
                vertex_reference[local_vertex][1]
            );
            const std::size_t node = triangle[local_vertex];
            accumulated[node].x += area * gradient.x;
            accumulated[node].y += area * gradient.y;
            accumulated_weights[node] += area;
        }
    }

    std::vector<Gradient2D> recovered(mesh.nodes.size());
    for (std::size_t node = 0; node < mesh.nodes.size(); ++node) {
        if (accumulated_weights[node] > 0.0) {
            recovered[node].x = accumulated[node].x / accumulated_weights[node];
            recovered[node].y = accumulated[node].y / accumulated_weights[node];
        }
        if (!std::isfinite(recovered[node].x) ||
            !std::isfinite(recovered[node].y)) {
            throw std::overflow_error("Recovered nodal gradient is non-finite.");
        }
    }

    H1GradientRecoveryIndicators2D result;
    result.element_errors.resize(mesh.elements.size(), 0.0);
    double total_squared = 0.0;

    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const Triangle2D& triangle = mesh.elements[element];
        const Point2D& p0 = mesh.nodes[triangle[0]];
        const Point2D& p1 = mesh.nodes[triangle[1]];
        const Point2D& p2 = mesh.nodes[triangle[2]];
        const double j00 = p1.x - p0.x;
        const double j01 = p2.x - p0.x;
        const double j10 = p1.y - p0.y;
        const double j11 = p2.y - p0.y;
        const double determinant = j00 * j11 - j01 * j10;
        const double absolute_determinant = std::abs(determinant);

        double indicator_squared = 0.0;
        for (const TriangleQuadraturePoint& quadrature : triangle_degree_four_quadrature) {
            const double x = p0.x + j00 * quadrature.r + j01 * quadrature.s;
            const double y = p0.y + j10 * quadrature.r + j11 * quadrature.s;
            const double k = conductivity ? conductivity(x, y) : 1.0;
            if (!std::isfinite(k) || k <= 0.0) {
                throw std::invalid_argument(
                    "Gradient recovery conductivity must be finite and positive."
                );
            }

            const Gradient2D finite_element_gradient = evaluate_element_gradient(
                mesh, dofs, order, dof_values, element,
                quadrature.r, quadrature.s
            );
            const auto shape = triangle_p1_shape_values(
                quadrature.r, quadrature.s
            );
            Gradient2D recovered_gradient;
            for (std::size_t local_vertex = 0; local_vertex < 3U; ++local_vertex) {
                const Gradient2D& nodal = recovered[triangle[local_vertex]];
                recovered_gradient.x += shape[local_vertex] * nodal.x;
                recovered_gradient.y += shape[local_vertex] * nodal.y;
            }

            const double difference_x =
                recovered_gradient.x - finite_element_gradient.x;
            const double difference_y =
                recovered_gradient.y - finite_element_gradient.y;
            indicator_squared += k *
                (difference_x * difference_x + difference_y * difference_y) *
                absolute_determinant * quadrature.weight;
        }

        if (!std::isfinite(indicator_squared) || indicator_squared < 0.0) {
            throw std::overflow_error("Gradient recovery indicator overflowed.");
        }
        result.element_errors[element] = std::sqrt(indicator_squared);
        total_squared += indicator_squared;
        if (!std::isfinite(total_squared)) {
            throw std::overflow_error("Global gradient recovery indicator overflowed.");
        }
    }

    result.total_error = std::sqrt(total_squared);
    return result;
}

H1GradientRecoveryIndicators2D
estimate_transient_heat_gradient_recovery_2d(
    const TransientHeat2DProblem& problem,
    const H1TriangleOrder order,
    const std::span<const double> temperature
) {
    if (!problem.conductivity) {
        throw std::invalid_argument(
            "Transient heat gradient recovery requires a conductivity field."
        );
    }
    return estimate_h1_gradient_recovery_error_2d(
        problem.mesh, order, temperature, problem.conductivity
    );
}

std::vector<std::size_t> mark_elements_dorfler_2d(
    const std::span<const double> element_errors,
    const double fraction
) {
    if (!std::isfinite(fraction) || fraction <= 0.0 || fraction > 1.0) {
        throw std::invalid_argument("Dorfler marking fraction must be in (0, 1].");
    }

    std::vector<std::pair<double, std::size_t>> ranked;
    ranked.reserve(element_errors.size());
    double total_squared = 0.0;
    for (std::size_t element = 0; element < element_errors.size(); ++element) {
        const double error = element_errors[element];
        if (!std::isfinite(error) || error < 0.0) {
            throw std::invalid_argument(
                "Dorfler marking indicators must be finite and nonnegative."
            );
        }
        const double squared = error * error;
        if (!std::isfinite(squared)) {
            throw std::overflow_error("Dorfler marking indicator square overflowed.");
        }
        total_squared += squared;
        if (!std::isfinite(total_squared)) {
            throw std::overflow_error("Dorfler marking total overflowed.");
        }
        if (squared > 0.0) {
            ranked.emplace_back(squared, element);
        }
    }

    if (total_squared == 0.0) {
        return {};
    }

    std::sort(
        ranked.begin(), ranked.end(),
        [](const auto& left, const auto& right) {
            if (left.first != right.first) {
                return left.first > right.first;
            }
            return left.second < right.second;
        }
    );

    const double target = fraction * total_squared;
    double accumulated = 0.0;
    std::vector<std::size_t> marked;
    for (const auto& [squared, element] : ranked) {
        marked.push_back(element);
        accumulated += squared;
        if (accumulated >= target) {
            break;
        }
    }
    std::sort(marked.begin(), marked.end());
    return marked;
}

} // namespace nova::fem
