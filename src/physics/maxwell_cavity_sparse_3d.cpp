#include "nova/physics/maxwell_cavity_sparse_3d.hpp"
#include "nova/linalg/conjugate_gradient.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::physics {
namespace {

using Edge = fem::Edge3D;
using linalg::CsrMatrix;
using linalg::Triplet;

struct GradientData {
    CsrMatrix incidence;
    CsrMatrix gram;
    std::size_t interior_nodes{};
};

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Maxwell cavity dot-product dimension mismatch.");
    }
    double sum = 0.0;
    double correction = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double product = a[i] * b[i];
        const double adjusted = product - correction;
        const double next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    if (!std::isfinite(sum)) {
        throw std::overflow_error("Maxwell cavity dot product overflow.");
    }
    return sum;
}

double norm(const std::vector<double>& x) {
    double value = 0.0;
    for (const double entry : x) {
        if (!std::isfinite(entry)) {
            throw std::overflow_error("Maxwell cavity vector is nonfinite.");
        }
        value = std::hypot(value, entry);
    }
    return value;
}

void add_scaled(std::vector<double>& y, const double scale,
                const std::vector<double>& x) {
    if (y.size() != x.size() || !std::isfinite(scale)) {
        throw std::invalid_argument("Invalid Maxwell cavity vector update.");
    }
    for (std::size_t i = 0; i < y.size(); ++i) {
        y[i] = std::fma(scale, x[i], y[i]);
        if (!std::isfinite(y[i])) {
            throw std::overflow_error("Maxwell cavity vector update overflow.");
        }
    }
}

double mass_norm(const std::vector<double>& x, const std::vector<double>& mx) {
    const double squared = dot(x, mx);
    if (squared < 0.0) {
        throw std::runtime_error("Maxwell cavity M norm became negative.");
    }
    return std::sqrt(squared);
}

double frobenius_norm(const CsrMatrix& matrix) {
    double value = 0.0;
    matrix.for_each_nonzero([&](std::size_t, std::size_t, const double entry) {
        value = std::hypot(value, entry);
    });
    return value;
}

void validate_symmetric(const CsrMatrix& matrix, const char* label) {
    if (matrix.rows() == 0U || matrix.rows() != matrix.columns()) {
        throw std::invalid_argument(std::string(label) + " must be non-empty and square.");
    }
    double scale = 0.0;
    matrix.for_each_nonzero([&](std::size_t, std::size_t, const double value) {
        scale = std::max(scale, std::abs(value));
    });
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * scale;
    bool symmetric = true;
    matrix.for_each_nonzero([&](const std::size_t row, const std::size_t column,
                                const double value) {
        if (std::abs(value - matrix.at(column, row)) > tolerance) {
            symmetric = false;
        }
    });
    if (!symmetric) {
        throw std::invalid_argument(std::string(label) + " must be symmetric.");
    }
}

std::vector<Edge> all_mesh_edges(const fem::TetrahedralMesh3D& mesh) {
    std::set<Edge> unique;
    for (const auto& cell : mesh.elements) {
        for (const auto& local : fem::tetrahedron_local_edges) {
            unique.insert({
                std::min(cell[local[0]], cell[local[1]]),
                std::max(cell[local[0]], cell[local[1]])
            });
        }
    }
    return std::vector<Edge>(unique.begin(), unique.end());
}

std::vector<std::size_t> all_free_edges(
    const fem::TetrahedralMesh3D& mesh, const std::vector<Edge>& edges) {
    std::set<Edge> pec;
    for (const auto& face : fem::extract_boundary_faces_3d(mesh)) {
        for (std::size_t i = 0; i < 3U; ++i) {
            for (std::size_t j = i + 1U; j < 3U; ++j) {
                pec.insert({
                    std::min(face[i], face[j]),
                    std::max(face[i], face[j])
                });
            }
        }
    }
    std::vector<std::size_t> free;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (!pec.contains(edges[i])) free.push_back(i);
    }
    return free;
}

GradientData build_gradient(
    const fem::TetrahedralMesh3D& mesh,
    const MaxwellCavityMatrices3D& matrices) {
    const std::size_t invalid = std::numeric_limits<std::size_t>::max();
    std::vector<unsigned char> boundary(mesh.nodes.size(), 0U);
    for (const auto& face : fem::extract_boundary_faces_3d(mesh)) {
        for (const auto node : face) boundary[node] = 1U;
    }

    std::vector<std::size_t> interior_index(mesh.nodes.size(), invalid);
    std::size_t interior_count = 0U;
    for (std::size_t node = 0; node < mesh.nodes.size(); ++node) {
        if (boundary[node] == 0U) interior_index[node] = interior_count++;
    }

    std::vector<Triplet> entries;
    if (matrices.free_edges.size() > std::numeric_limits<std::size_t>::max() / 2U) {
        throw std::length_error("Maxwell cavity gradient incidence is too large.");
    }
    entries.reserve(2U * matrices.free_edges.size());
    for (std::size_t row = 0; row < matrices.free_edges.size(); ++row) {
        const auto& edge = matrices.edges[matrices.free_edges[row]];
        if (interior_index[edge[0]] != invalid) {
            entries.push_back({row, interior_index[edge[0]], -1.0});
        }
        if (interior_index[edge[1]] != invalid) {
            entries.push_back({row, interior_index[edge[1]], 1.0});
        }
    }
    CsrMatrix incidence = CsrMatrix::from_triplets(
        matrices.free_edges.size(), interior_count, entries);
    if (interior_count == 0U) {
        return {std::move(incidence), CsrMatrix{}, 0U};
    }

    struct NodeSign { std::size_t node{}; double sign{}; };
    std::vector<std::vector<NodeSign>> edge_nodes(matrices.free_edges.size());
    incidence.for_each_nonzero([&](const std::size_t edge,
                                   const std::size_t node,
                                   const double sign) {
        edge_nodes[edge].push_back({node, sign});
    });

    if (matrices.mass.non_zeros() > std::numeric_limits<std::size_t>::max() / 4U) {
        throw std::length_error("Maxwell cavity gradient Gram matrix is too large.");
    }
    std::vector<Triplet> gram_entries;
    gram_entries.reserve(4U * matrices.mass.non_zeros());
    matrices.mass.for_each_nonzero([&](const std::size_t row,
                                       const std::size_t column,
                                       const double value) {
        for (const auto& left : edge_nodes[row]) {
            for (const auto& right : edge_nodes[column]) {
                const double entry = left.sign * value * right.sign;
                if (!std::isfinite(entry)) {
                    throw std::overflow_error("Maxwell cavity gradient Gram entry overflow.");
                }
                if (entry != 0.0) gram_entries.push_back({left.node, right.node, entry});
            }
        }
    });
    CsrMatrix gram = CsrMatrix::from_triplets(
        interior_count, interior_count, gram_entries);
    for (const double value : gram.diagonal()) {
        if (!(value > 0.0) || !std::isfinite(value)) {
            throw std::invalid_argument("Discrete H1_0 gradient Gram matrix is not positive definite.");
        }
    }
    return {std::move(incidence), std::move(gram), interior_count};
}

std::uint64_t splitmix64(std::uint64_t& state) {
    std::uint64_t value = (state += 0x9e3779b97f4a7c15ULL);
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

void project_gradient(
    std::vector<double>& vector,
    const GradientData& gradient,
    const CsrMatrix& mass,
    const linalg::SolverOptions& options,
    MaxwellCavitySparseResult3D& result) {
    if (gradient.interior_nodes == 0U) return;

    const auto mv = mass.multiply(vector);
    std::vector<double> rhs(gradient.interior_nodes, 0.0);
    gradient.incidence.for_each_nonzero(
        [&](const std::size_t edge, const std::size_t node, const double sign) {
            rhs[node] += sign * mv[edge];
        });
    if (norm(rhs) == 0.0) return;

    const auto solve = linalg::solve_conjugate_gradient(
        gradient.gram, rhs, options);
    result.projection_solve_iterations += solve.iterations;
    if (!solve.converged) {
        throw std::runtime_error(
            "Maxwell discrete-gradient projection failed: " + solve.message);
    }
    const auto correction = gradient.incidence.multiply(solve.solution);
    for (std::size_t i = 0; i < vector.size(); ++i) {
        vector[i] -= correction[i];
        if (!std::isfinite(vector[i])) {
            throw std::overflow_error("Maxwell gradient projection overflow.");
        }
    }
}

} // namespace

MaxwellCavitySparseResult3D solve_maxwell_cavity_shift_invert_3d(
    const fem::TetrahedralMesh3D& mesh,
    const MaxwellCavityMatrices3D& matrices,
    const MaxwellCavitySparseOptions3D& options) {
    if (options.max_modes == 0U || options.max_subspace == 0U ||
        options.max_subspace > 512U || options.max_shift_iterations == 0U ||
        options.max_projection_iterations == 0U ||
        !std::isfinite(options.shift_scale) || !(options.shift_scale > 0.0) ||
        !std::isfinite(options.relative_tolerance) ||
        !(options.relative_tolerance > 0.0) || options.relative_tolerance >= 1.0 ||
        !std::isfinite(options.shift_solve_relative_tolerance) ||
        !(options.shift_solve_relative_tolerance > 0.0) ||
        options.shift_solve_relative_tolerance >= 1.0 ||
        !std::isfinite(options.projection_relative_tolerance) ||
        !(options.projection_relative_tolerance > 0.0) ||
        options.projection_relative_tolerance >= 1.0 ||
        !std::isfinite(options.nullspace_relative_tolerance) ||
        !(options.nullspace_relative_tolerance > 0.0) ||
        options.nullspace_relative_tolerance >= 1.0 ||
        !std::isfinite(options.breakdown_tolerance) ||
        !(options.breakdown_tolerance > 0.0) || options.breakdown_tolerance >= 1.0) {
        throw std::invalid_argument("Invalid sparse Maxwell cavity options.");
    }

    mesh.validate();
    const auto expected_edges = all_mesh_edges(mesh);
    const auto expected_free = all_free_edges(mesh, expected_edges);
    if (matrices.edges != expected_edges || matrices.free_edges != expected_free ||
        matrices.free_edges.size() != matrices.stiffness.rows() ||
        matrices.stiffness.rows() != matrices.mass.rows() ||
        matrices.stiffness.columns() != matrices.mass.columns()) {
        throw std::invalid_argument("Maxwell cavity matrices do not match the PEC mesh topology.");
    }
    validate_symmetric(matrices.stiffness, "Maxwell cavity stiffness");
    validate_symmetric(matrices.mass, "Maxwell cavity mass");

    const auto kd = matrices.stiffness.diagonal();
    const auto md = matrices.mass.diagonal();
    double spectral_scale = 0.0;
    for (std::size_t i = 0; i < kd.size(); ++i) {
        if (!(md[i] > 0.0) || !std::isfinite(md[i]) ||
            !std::isfinite(kd[i]) || kd[i] < 0.0) {
            throw std::invalid_argument("Maxwell cavity requires non-negative K and positive M diagonals.");
        }
        spectral_scale = std::max(spectral_scale, kd[i] / md[i]);
    }
    if (!(spectral_scale > 0.0) || !std::isfinite(spectral_scale)) {
        throw std::invalid_argument("Maxwell cavity has no positive stiffness scale.");
    }

    MaxwellCavitySparseResult3D result;
    const GradientData gradient = build_gradient(mesh, matrices);
    result.gradient_nullspace_dimension = gradient.interior_nodes;
    const std::size_t n = matrices.stiffness.rows();
    if (gradient.interior_nodes >= n || n - gradient.interior_nodes < options.max_modes) {
        throw std::invalid_argument("Too few non-gradient Maxwell edge degrees of freedom.");
    }
    const std::size_t complement_dimension = n - gradient.interior_nodes;
    const std::size_t subspace_limit = std::min(options.max_subspace, complement_dimension);
    if (subspace_limit < options.max_modes) {
        throw std::invalid_argument("Krylov subspace is smaller than the requested mode count.");
    }
    const std::size_t start_count = std::min(subspace_limit, options.max_modes + 2U);

    result.shift = options.shift_scale * spectral_scale;
    result.nullspace_cutoff = options.nullspace_relative_tolerance * spectral_scale;
    if (!(result.shift > 0.0) || !std::isfinite(result.shift) ||
        !(result.nullspace_cutoff > 0.0) || !std::isfinite(result.nullspace_cutoff)) {
        throw std::overflow_error("Maxwell cavity shift or nullspace cutoff overflow.");
    }
    const CsrMatrix shifted = CsrMatrix::linear_combination(
        matrices.stiffness, 1.0, matrices.mass, result.shift);

    linalg::SolverOptions shift_options;
    shift_options.max_iterations = options.max_shift_iterations;
    shift_options.relative_tolerance = options.shift_solve_relative_tolerance;
    shift_options.absolute_tolerance = 0.0;
    shift_options.preconditioner = linalg::PreconditionerType::Jacobi;

    linalg::SolverOptions projection_options;
    projection_options.max_iterations = options.max_projection_iterations;
    projection_options.relative_tolerance = options.projection_relative_tolerance;
    projection_options.absolute_tolerance = 0.0;
    projection_options.preconditioner = linalg::PreconditionerType::Jacobi;

    std::vector<std::vector<double>> basis;
    std::vector<std::vector<double>> mass_basis;
    std::vector<std::vector<double>> operator_mass_basis;
    basis.reserve(subspace_limit);
    mass_basis.reserve(subspace_limit);
    operator_mass_basis.reserve(subspace_limit);

    std::uint64_t random_state = 0x4e4f56415f4d4158ULL;
    for (std::size_t seed = 0; seed < start_count; ++seed) {
        bool accepted = false;
        for (std::size_t attempt = 0; attempt < 16U * start_count + 8U; ++attempt) {
            std::vector<double> candidate(n);
            for (double& value : candidate) {
                value = (splitmix64(random_state) & 1ULL) != 0ULL ? 1.0 : -1.0;
            }
            project_gradient(candidate, gradient, matrices.mass, projection_options, result);
            auto mc = matrices.mass.multiply(candidate);
            const double initial_norm = mass_norm(candidate, mc);
            if (!(initial_norm > 0.0)) continue;

            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t j = 0; j < basis.size(); ++j) {
                    const double coefficient = dot(basis[j], mc);
                    add_scaled(candidate, -coefficient, basis[j]);
                    add_scaled(mc, -coefficient, mass_basis[j]);
                }
            }
            mc = matrices.mass.multiply(candidate);
            const double residual_norm = mass_norm(candidate, mc);
            if (!(residual_norm > options.breakdown_tolerance * initial_norm)) continue;

            const double inverse_norm = 1.0 / residual_norm;
            for (double& value : candidate) value *= inverse_norm;
            mc = matrices.mass.multiply(candidate);
            basis.push_back(std::move(candidate));
            mass_basis.push_back(std::move(mc));
            accepted = true;
            break;
        }
        if (!accepted) {
            throw std::runtime_error("Unable to construct independent M-orthonormal Maxwell starting vectors.");
        }
    }

    for (std::size_t current = 0; current < basis.size() &&
         current < subspace_limit; ++current) {
        const auto rhs = matrices.mass.multiply(basis[current]);
        const auto solve = linalg::solve_conjugate_gradient(shifted, rhs, shift_options);
        result.shift_solve_iterations += solve.iterations;
        if (!solve.converged) {
            throw std::runtime_error("Maxwell shift-invert solve failed: " + solve.message);
        }

        auto transformed = solve.solution;
        project_gradient(transformed, gradient, matrices.mass, projection_options, result);
        auto mt = matrices.mass.multiply(transformed);
        const double initial_norm = mass_norm(transformed, mt);
        if (!(initial_norm > 0.0)) {
            throw std::runtime_error("Maxwell shift-invert produced a zero vector.");
        }

        auto candidate = transformed;
        auto mc = mt;
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t j = 0; j < basis.size(); ++j) {
                const double coefficient = dot(basis[j], mc);
                add_scaled(candidate, -coefficient, basis[j]);
                add_scaled(mc, -coefficient, mass_basis[j]);
            }
        }
        mc = matrices.mass.multiply(candidate);
        const double residual_norm = mass_norm(candidate, mc);
        if (residual_norm > options.breakdown_tolerance * initial_norm &&
            basis.size() < subspace_limit) {
            const double inverse_norm = 1.0 / residual_norm;
            for (double& value : candidate) value *= inverse_norm;
            mc = matrices.mass.multiply(candidate);
            basis.push_back(std::move(candidate));
            mass_basis.push_back(std::move(mc));
        }
        operator_mass_basis.push_back(std::move(mt));
    }

    if (operator_mass_basis.size() != basis.size()) {
        throw std::runtime_error("Maxwell Krylov basis was not fully processed.");
    }
    result.subspace_dimension = basis.size();
    result.operator_applications = operator_mass_basis.size();

    std::vector<Triplet> projected_entries;
    if (basis.size() > std::numeric_limits<std::size_t>::max() / basis.size()) {
        throw std::length_error("Maxwell projected eigenproblem is too large.");
    }
    projected_entries.reserve(basis.size() * basis.size());
    for (std::size_t row = 0; row < basis.size(); ++row) {
        for (std::size_t column = row; column < basis.size(); ++column) {
            const double forward = dot(basis[row], operator_mass_basis[column]);
            const double reverse = dot(basis[column], operator_mass_basis[row]);
            const double value = 0.5 * forward + 0.5 * reverse;
            if (!std::isfinite(value)) {
                throw std::overflow_error("Maxwell projected operator is nonfinite.");
            }
            if (value != 0.0) {
                projected_entries.push_back({row, column, value});
                if (row != column) projected_entries.push_back({column, row, value});
            }
        }
    }

    const CsrMatrix projected_operator = CsrMatrix::from_triplets(
        basis.size(), basis.size(), projected_entries);
    std::vector<Triplet> identity_entries;
    identity_entries.reserve(basis.size());
    for (std::size_t i = 0; i < basis.size(); ++i) identity_entries.push_back({i, i, 1.0});
    const CsrMatrix identity = CsrMatrix::from_triplets(
        basis.size(), basis.size(), identity_entries);
    linalg::SymmetricEigenOptions dense_options;
    dense_options.max_dimension = 512U;
    dense_options.max_sweeps = 120U;
    dense_options.relative_tolerance = 1.0e-13;
    const auto projected = linalg::solve_symmetric_generalized_eigen(
        projected_operator, identity, dense_options);

    const double k_norm = frobenius_norm(matrices.stiffness);
    const double m_norm = frobenius_norm(matrices.mass);
    if (!std::isfinite(k_norm) || !std::isfinite(m_norm) || !(m_norm > 0.0)) {
        throw std::overflow_error("Maxwell cavity matrix norm is outside the supported range.");
    }
    bool residuals_converged = projected.converged;
    for (std::size_t offset = 0;
         offset < projected.eigenpairs.size() && result.modes.size() < options.max_modes;
         ++offset) {
        const auto& ritz = projected.eigenpairs[projected.eigenpairs.size() - 1U - offset];
        const double mu = ritz.eigenvalue;
        if (!(mu > 0.0) || !std::isfinite(mu)) {
            throw std::runtime_error("Maxwell shift-invert Ritz value is not positive.");
        }
        const double lambda = 1.0 / mu - result.shift;
        if (!std::isfinite(lambda)) {
            throw std::overflow_error("Maxwell eigenvalue conversion overflow.");
        }
        if (lambda < -result.nullspace_cutoff) {
            throw std::runtime_error("Maxwell curl operator has a negative physical eigenvalue.");
        }
        if (lambda <= result.nullspace_cutoff) {
            ++result.additional_zero_modes_detected;
            continue;
        }

        std::vector<double> x(n, 0.0);
        for (std::size_t j = 0; j < basis.size(); ++j) {
            add_scaled(x, ritz.eigenvector[j], basis[j]);
        }
        auto mx = matrices.mass.multiply(x);
        const double x_norm = mass_norm(x, mx);
        if (!(x_norm > 0.0) || !std::isfinite(x_norm)) {
            throw std::runtime_error("Maxwell Ritz vector has invalid M norm.");
        }
        const double inverse_norm = 1.0 / x_norm;
        for (double& value : x) value *= inverse_norm;

        const auto kx = matrices.stiffness.multiply(x);
        mx = matrices.mass.multiply(x);
        std::vector<double> residual(n);
        for (std::size_t i = 0; i < n; ++i) residual[i] = kx[i] - lambda * mx[i];
        const double numerator = norm(residual);
        const double denominator = (k_norm + std::abs(lambda) * m_norm) * norm(x);
        if (!std::isfinite(denominator)) {
            throw std::overflow_error("Maxwell cavity residual normalization overflow.");
        }
        const double relative_residual = denominator > 0.0
            ? numerator / denominator : numerator;
        if (!std::isfinite(relative_residual)) {
            throw std::overflow_error("Maxwell Ritz residual is nonfinite.");
        }
        if (relative_residual > options.relative_tolerance) residuals_converged = false;

        const auto largest = std::max_element(
            x.begin(), x.end(), [](const double a, const double b) {
                return std::abs(a) < std::abs(b);
            });
        if (largest != x.end() && *largest < 0.0) {
            for (double& value : x) value = -value;
        }

        MaxwellCavityMode3D mode;
        mode.vacuum_wavenumber_squared = lambda;
        mode.frequency_hz = 299792458.0 /
            (2.0 * std::numbers::pi) * std::sqrt(lambda);
        if (!std::isfinite(mode.frequency_hz)) {
            throw std::overflow_error("Maxwell cavity frequency overflow.");
        }
        mode.relative_residual = relative_residual;
        mode.edge_field.resize(matrices.edges.size(), 0.0);
        for (std::size_t i = 0; i < matrices.free_edges.size(); ++i) {
            mode.edge_field[matrices.free_edges[i]] = x[i];
        }
        result.modes.push_back(std::move(mode));
    }

    result.nullspace_dimension =
        result.gradient_nullspace_dimension + result.additional_zero_modes_detected;
    result.converged = projected.converged &&
        result.modes.size() == options.max_modes && residuals_converged;
    return result;
}

} // namespace nova::physics
