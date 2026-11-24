#include "nova/fem/residual_error_estimator_3d.hpp"

#include "nova/fem/tetrahedron_h1.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace nova::fem {
namespace {

struct GaussPoint1D {
    double coordinate;
    double weight;
};

// Five-point Gauss-Legendre quadrature on [0,1], degree-nine exact.
constexpr std::array<GaussPoint1D, 5> gauss_legendre_5{{
    {0.046910077030668, 0.118463442528095},
    {0.230765344947158, 0.239314335249683},
    {0.5,               0.284444444444444},
    {0.769234655052842, 0.239314335249683},
    {0.953089922969332, 0.118463442528095}
}};

[[nodiscard]] double checked_norm(
    const long double square,
    const char* message
) {
    if (!(square >= 0.0L) || !std::isfinite(square)) {
        throw std::overflow_error(message);
    }
    const long double magnitude = std::sqrt(square);
    if (!std::isfinite(magnitude) ||
        magnitude > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error(message);
    }
    return static_cast<double>(magnitude);
}

[[nodiscard]] double tetrahedron_diameter(
    const TetrahedralMesh3D& mesh,
    const Tetrahedron3D& cell
) {
    double diameter = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = i + 1; j < 4; ++j) {
            const auto& a = mesh.nodes[cell[i]];
            const auto& b = mesh.nodes[cell[j]];
            const double edge = std::hypot(
                b.x - a.x, b.y - a.y, b.z - a.z);
            if (!std::isfinite(edge) || !(edge > 0.0)) {
                throw std::invalid_argument(
                    "Residual estimator encountered invalid tetrahedron geometry.");
            }
            diameter = std::max(diameter, edge);
        }
    }
    if (!std::isfinite(diameter) || !(diameter > 0.0)) {
        throw std::invalid_argument(
            "Residual estimator encountered invalid tetrahedron diameter.");
    }
    return diameter;
}

[[nodiscard]] Point3D map_reference_tetrahedron_point(
    const TetrahedralMesh3D& mesh,
    const Tetrahedron3D& cell,
    const double r,
    const double s,
    const double t
) {
    const auto& p0 = mesh.nodes[cell[0]];
    const auto& p1 = mesh.nodes[cell[1]];
    const auto& p2 = mesh.nodes[cell[2]];
    const auto& p3 = mesh.nodes[cell[3]];
    return {
        p0.x + r * (p1.x - p0.x) + s * (p2.x - p0.x) + t * (p3.x - p0.x),
        p0.y + r * (p1.y - p0.y) + s * (p2.y - p0.y) + t * (p3.y - p0.y),
        p0.z + r * (p1.z - p0.z) + s * (p2.z - p0.z) + t * (p3.z - p0.z)
    };
}

[[nodiscard]] long double integrate_source_square(
    const TetrahedralMesh3D& mesh,
    const Tetrahedron3D& cell,
    const double jacobian,
    const std::function<double(double, double, double)>& source
) {
    long double integral = 0.0L;
    for (const auto& q_u : gauss_legendre_5) {
        const double u = q_u.coordinate;
        const double one_minus_u = 1.0 - u;
        for (const auto& q_v : gauss_legendre_5) {
            const double v = q_v.coordinate;
            const double one_minus_v = 1.0 - v;
            for (const auto& q_w : gauss_legendre_5) {
                // Duffy map [0,1]^3 to the reference tetrahedron:
                // r=u, s=(1-u)v, t=(1-u)(1-v)w.
                const double r = u;
                const double s = one_minus_u * v;
                const double t = one_minus_u * one_minus_v * q_w.coordinate;
                const Point3D point =
                    map_reference_tetrahedron_point(mesh, cell, r, s, t);
                const double value = source(point.x, point.y, point.z);
                if (!std::isfinite(value)) {
                    throw std::invalid_argument(
                        "Residual estimator source returned a non-finite value.");
                }
                const long double duffy_jacobian =
                    static_cast<long double>(one_minus_u) * one_minus_u *
                    one_minus_v;
                const long double weight =
                    static_cast<long double>(jacobian) *
                    q_u.weight * q_v.weight * q_w.weight * duffy_jacobian;
                const long double f = value;
                integral += weight * f * f;
            }
        }
    }
    if (!(integral >= 0.0L) || !std::isfinite(integral)) {
        throw std::overflow_error(
            "Residual estimator source integral is not representable.");
    }
    return integral;
}


struct ResidualFaceOwner3D {
    std::size_t cell{};
    std::size_t opposite{};
    std::size_t occurrences{};
};

struct LongVector3D {
    long double x{}, y{}, z{};
};

[[nodiscard]] LongVector3D vector_between(const Point3D& a, const Point3D& b) {
    return {static_cast<long double>(b.x) - a.x,
            static_cast<long double>(b.y) - a.y,
            static_cast<long double>(b.z) - a.z};
}

[[nodiscard]] LongVector3D vector_cross(
    const LongVector3D& a, const LongVector3D& b
) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

[[nodiscard]] long double vector_dot(
    const LongVector3D& a, const LongVector3D& b
) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

[[nodiscard]] double checked_boundary_value(const long double value) {
    if (!std::isfinite(value) ||
        std::abs(value) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("3D natural-boundary residual is not representable.");
    }
    return static_cast<double>(value);
}

struct FaceQuadraturePoint3D {
    std::array<double, 3> barycentric;
    double normalized_weight;
};

// Positive degree-four triangle quadrature, weights normalized to sum to 1.
constexpr std::array<FaceQuadraturePoint3D, 6> face_quadrature{{
    {{{0.445948490915965, 0.445948490915965, 0.108103018168070}}, 0.223381589678010},
    {{{0.445948490915965, 0.108103018168070, 0.445948490915965}}, 0.223381589678010},
    {{{0.108103018168070, 0.445948490915965, 0.445948490915965}}, 0.223381589678010},
    {{{0.091576213509771, 0.091576213509771, 0.816847572980459}}, 0.109951743655322},
    {{{0.091576213509771, 0.816847572980459, 0.091576213509771}}, 0.109951743655322},
    {{{0.816847572980459, 0.091576213509771, 0.091576213509771}}, 0.109951743655322}
}};

struct SelectedNaturalFace3D {
    BoundaryFace3D face{};
    const NeumannResidualFace3D* neumann{};
    const RobinResidualFace3D* robin{};
};

[[nodiscard]] std::vector<long double> integrate_natural_residuals_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const P1TetrahedralFlux3D> fields,
    std::span<const NeumannResidualFace3D> natural_faces,
    std::span<const RobinResidualFace3D> robin_faces,
    std::span<const double> nodal_values,
    std::vector<BoundaryFluxResidualFace3D>& report
) {
    std::vector<long double> squared(mesh.elements.size(), 0.0L);
    if (natural_faces.empty() && robin_faces.empty()) return squared;

    std::vector<SelectedNaturalFace3D> selected_faces;
    selected_faces.reserve(natural_faces.size() + robin_faces.size());
    for (const auto& face : natural_faces) {
        selected_faces.push_back({face.face, &face, nullptr});
    }
    for (const auto& face : robin_faces) {
        selected_faces.push_back({face.face, nullptr, &face});
    }

    std::map<BoundaryFace3D, ResidualFaceOwner3D> owners;
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const auto& tet = mesh.elements[element];
        for (std::size_t omitted = 0; omitted < 4; ++omitted) {
            BoundaryFace3D face{};
            std::size_t position = 0;
            for (std::size_t vertex = 0; vertex < 4; ++vertex) {
                if (vertex != omitted) face[position++] = tet[vertex];
            }
            std::sort(face.begin(), face.end());
            const auto [it, inserted] = owners.try_emplace(
                face, ResidualFaceOwner3D{element, tet[omitted], 1U});
            if (!inserted) ++it->second.occurrences;
        }
    }

    std::map<BoundaryFace3D, bool> visited;
    report.reserve(selected_faces.size());
    for (const auto& selected : selected_faces) {
        if (selected.neumann && !selected.neumann->outward_flux) {
            throw std::invalid_argument(
                "Natural-boundary residual requires an outward-flux callback.");
        }
        if (selected.robin &&
            (!selected.robin->transfer_coefficient || !selected.robin->ambient_value)) {
            throw std::invalid_argument(
                "Robin residual requires transfer and ambient callbacks.");
        }
        BoundaryFace3D face = selected.face;
        for (const auto vertex : face) {
            if (vertex >= mesh.nodes.size()) {
                throw std::out_of_range(
                    "Natural-boundary residual references a missing vertex.");
            }
        }
        std::sort(face.begin(), face.end());
        const auto found = owners.find(face);
        if (found == owners.end() || found->second.occurrences != 1U) {
            throw std::invalid_argument(
                "Natural-boundary residual face must be an exterior tetrahedral face.");
        }
        if (!visited.emplace(face, true).second) {
            throw std::invalid_argument("Natural-boundary residual face is duplicated.");
        }

        const auto& a = mesh.nodes[face[0]];
        const auto& b = mesh.nodes[face[1]];
        const auto& c = mesh.nodes[face[2]];
        const auto& opposite = mesh.nodes[found->second.opposite];
        const auto ab = vector_between(a, b);
        const auto ac = vector_between(a, c);
        auto normal_area_twice = vector_cross(ab, ac);
        const long double length_twice = std::hypot(
            normal_area_twice.x, normal_area_twice.y, normal_area_twice.z);
        const long double area = 0.5L * length_twice;
        const long double diameter = std::max({
            std::hypot(ab.x, ab.y, ab.z),
            std::hypot(ac.x, ac.y, ac.z),
            std::hypot(static_cast<long double>(c.x) - b.x,
                       static_cast<long double>(c.y) - b.y,
                       static_cast<long double>(c.z) - b.z)
        });
        if (!(area > 0.0L) || !(diameter > 0.0L) ||
            !std::isfinite(area) || !std::isfinite(diameter)) {
            throw std::overflow_error(
                "Natural-boundary residual face has invalid geometry.");
        }
        if (vector_dot(normal_area_twice, vector_between(a, opposite)) > 0.0L) {
            normal_area_twice.x = -normal_area_twice.x;
            normal_area_twice.y = -normal_area_twice.y;
            normal_area_twice.z = -normal_area_twice.z;
        }
        const auto& flux = fields[found->second.cell].flux;
        const LongVector3D q{flux.x, flux.y, flux.z};
        const long double actual_outward =
            vector_dot(q, normal_area_twice) / length_twice;
        long double weighted_mismatch_square = 0.0L;
        for (const auto& qp : face_quadrature) {
            const auto& weight = qp.barycentric;
            const long double x = weight[0] * a.x + weight[1] * b.x + weight[2] * c.x;
            const long double y = weight[0] * a.y + weight[1] * b.y + weight[2] * c.y;
            const long double z = weight[0] * a.z + weight[1] * b.z + weight[2] * c.z;
            const double px = checked_boundary_value(x);
            const double py = checked_boundary_value(y);
            const double pz = checked_boundary_value(z);
            long double prescribed = 0.0L;
            if (selected.neumann) {
                const double flux_value = selected.neumann->outward_flux(px, py, pz);
                if (!std::isfinite(flux_value)) {
                    throw std::invalid_argument(
                        "Natural-boundary flux callback returned a non-finite value.");
                }
                prescribed = flux_value;
            } else {
                const double h = selected.robin->transfer_coefficient(px, py, pz);
                const double ambient = selected.robin->ambient_value(px, py, pz);
                if (!std::isfinite(h) || h < 0.0 || !std::isfinite(ambient)) {
                    throw std::invalid_argument(
                        "Robin residual requires finite nonnegative h and finite ambient.");
                }
                const long double trace =
                    static_cast<long double>(weight[0]) * nodal_values[face[0]] +
                    static_cast<long double>(weight[1]) * nodal_values[face[1]] +
                    static_cast<long double>(weight[2]) * nodal_values[face[2]];
                prescribed = static_cast<long double>(h) * (trace - ambient);
                if (!std::isfinite(prescribed)) {
                    throw std::overflow_error("Robin residual prescribed flux overflow.");
                }
            }
            const long double mismatch = actual_outward - prescribed;
            weighted_mismatch_square += qp.normalized_weight * mismatch * mismatch;
        }
        const long double mismatch_l2_square = area * weighted_mismatch_square;
        const long double contribution = diameter * mismatch_l2_square;
        if (!(contribution >= 0.0L) || !std::isfinite(contribution)) {
            throw std::overflow_error("Natural-boundary residual exceeds finite range.");
        }
        squared[found->second.cell] += contribution;
        report.push_back({
            face, found->second.cell,
            checked_boundary_value(area),
            checked_boundary_value(diameter),
            checked_norm(mismatch_l2_square,
                "Natural-boundary L2 mismatch exceeds finite range."),
            selected.robin ? BoundaryResidualKind3D::Robin
                           : BoundaryResidualKind3D::Neumann
        });
    }
    return squared;
}

} // namespace

DiffusionResidualIndicators3D estimate_diffusion_residual_error_3d(
    const TetrahedralMesh3D& mesh,
    const std::span<const double> nodal_values,
    const std::span<const P1TetrahedralFlux3D> cell_fields,
    const std::function<double(double, double, double)>& source,
    const std::span<const NeumannResidualFace3D> natural_boundary_faces,
    const std::span<const RobinResidualFace3D> robin_boundary_faces,
    const std::span<const double> cell_source
) {
    mesh.validate();
    if (nodal_values.size() != mesh.nodes.size()) {
        throw std::invalid_argument(
            "3D residual estimator nodal field size must match the mesh node count.");
    }
    if (cell_fields.size() != mesh.elements.size()) {
        throw std::invalid_argument(
            "3D residual estimator flux count must match the element count.");
    }
    if (!source) {
        throw std::invalid_argument(
            "3D residual estimator requires a source function.");
    }
    if (!cell_source.empty() && cell_source.size() != mesh.elements.size()) {
        throw std::invalid_argument(
            "3D residual estimator cell-source count must match the element count.");
    }
    for (const double value : cell_source) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "3D residual estimator cell sources must be finite.");
        }
    }
    for (const double value : nodal_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "3D residual estimator nodal field must be finite.");
        }
    }
    for (const auto& field : cell_fields) {
        if (!std::isfinite(field.gradient.x) || !std::isfinite(field.gradient.y) ||
            !std::isfinite(field.gradient.z) || !std::isfinite(field.flux.x) ||
            !std::isfinite(field.flux.y) || !std::isfinite(field.flux.z)) {
            throw std::invalid_argument(
                "3D residual estimator requires finite reconstructed cell fields.");
        }
    }

    std::vector<long double> cell_residual_squared(mesh.elements.size(), 0.0L);
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const auto& cell = mesh.elements[element];
        const auto basis =
            evaluate_tetrahedron_p1_physical_gradients_3d(mesh, element);
        const double jacobian = std::abs(basis.signed_jacobian_determinant);
        if (!std::isfinite(jacobian) || !(jacobian > 0.0)) {
            throw std::invalid_argument(
                "3D residual estimator encountered invalid element volume.");
        }
        const double diameter = tetrahedron_diameter(mesh, cell);
        const double cell_source_value =
            cell_source.empty() ? 0.0 : cell_source[element];
        const std::function<double(double, double, double)> combined_source =
            [&source, cell_source_value](double x, double y, double z) {
                return source(x, y, z) + cell_source_value;
            };
        const long double source_integral =
            integrate_source_square(mesh, cell, jacobian, combined_source);
        const long double h = diameter;
        cell_residual_squared[element] = h * h * source_integral;
        if (!(cell_residual_squared[element] >= 0.0L) ||
            !std::isfinite(cell_residual_squared[element])) {
            throw std::overflow_error(
                "3D cell residual indicator is not representable.");
        }
    }

    DiffusionResidualIndicators3D result;
    result.flux_jumps = evaluate_interior_flux_jumps_3d(mesh, cell_fields);
    const auto boundary_squared = integrate_natural_residuals_3d(
        mesh, cell_fields, natural_boundary_faces, robin_boundary_faces,
        nodal_values, result.boundary_residual_faces);
    result.cell_residual_errors.resize(mesh.elements.size(), 0.0);
    result.boundary_residual_errors.resize(mesh.elements.size(), 0.0);
    result.element_errors.resize(mesh.elements.size(), 0.0);

    long double total_squared = 0.0L;
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const long double face_squared =
            static_cast<long double>(result.flux_jumps.cell_indicators[element]) *
            result.flux_jumps.cell_indicators[element];
        const long double combined_squared =
            cell_residual_squared[element] + face_squared + boundary_squared[element];
        result.cell_residual_errors[element] = checked_norm(
            cell_residual_squared[element],
            "3D cell residual indicator exceeds the finite range.");
        result.boundary_residual_errors[element] = checked_norm(
            boundary_squared[element],
            "3D boundary residual indicator exceeds the finite range.");
        result.element_errors[element] = checked_norm(
            combined_squared,
            "3D element residual indicator exceeds the finite range.");
        total_squared += combined_squared;
    }
    result.total_error = checked_norm(
        total_squared, "3D total residual indicator exceeds the finite range.");
    return result;
}

DiffusionResidualIndicators3D estimate_diffusion_residual_error_3d(
    const TetrahedralMesh3D& mesh,
    const std::span<const double> nodal_values,
    const std::span<const double> element_conductivity,
    const std::function<double(double, double, double)>& source,
    const std::span<const NeumannResidualFace3D> natural_boundary_faces,
    const std::span<const RobinResidualFace3D> robin_boundary_faces
) {
    mesh.validate();
    if (nodal_values.size() != mesh.nodes.size()) {
        throw std::invalid_argument(
            "3D residual estimator nodal field size must match the mesh node count.");
    }
    if (element_conductivity.size() != mesh.elements.size()) {
        throw std::invalid_argument(
            "3D residual estimator conductivity count must match the element count.");
    }
    if (!source) {
        throw std::invalid_argument(
            "3D residual estimator requires a source function.");
    }
    for (const double value : nodal_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "3D residual estimator nodal field must be finite.");
        }
    }
    for (const double conductivity : element_conductivity) {
        if (!std::isfinite(conductivity) || !(conductivity > 0.0)) {
            throw std::invalid_argument(
                "3D residual estimator conductivity must be finite and positive.");
        }
    }

    std::vector<P1TetrahedralFlux3D> cell_fields;
    cell_fields.reserve(mesh.elements.size());
    std::vector<long double> cell_residual_squared(mesh.elements.size(), 0.0L);

    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const auto& cell = mesh.elements[element];
        const auto basis =
            evaluate_tetrahedron_p1_physical_gradients_3d(mesh, element);
        const double jacobian = std::abs(basis.signed_jacobian_determinant);
        if (!std::isfinite(jacobian) || !(jacobian > 0.0)) {
            throw std::invalid_argument(
                "3D residual estimator encountered invalid element volume.");
        }

        Point3D gradient{};
        for (std::size_t local = 0; local < 4; ++local) {
            const double value = nodal_values[cell[local]];
            gradient.x += value * basis.gradients[local][0];
            gradient.y += value * basis.gradients[local][1];
            gradient.z += value * basis.gradients[local][2];
        }
        const double conductivity = element_conductivity[element];
        const Point3D flux{
            -conductivity * gradient.x,
            -conductivity * gradient.y,
            -conductivity * gradient.z
        };
        if (!std::isfinite(gradient.x) || !std::isfinite(gradient.y) ||
            !std::isfinite(gradient.z) || !std::isfinite(flux.x) ||
            !std::isfinite(flux.y) || !std::isfinite(flux.z)) {
            throw std::overflow_error(
                "3D residual estimator reconstructed a non-finite cell field.");
        }
        cell_fields.push_back({gradient, flux});

        const double diameter = tetrahedron_diameter(mesh, cell);
        const long double source_integral =
            integrate_source_square(mesh, cell, jacobian, source);
        const long double h = diameter;
        cell_residual_squared[element] = h * h * source_integral;
        if (!(cell_residual_squared[element] >= 0.0L) ||
            !std::isfinite(cell_residual_squared[element])) {
            throw std::overflow_error(
                "3D cell residual indicator is not representable.");
        }
    }

    DiffusionResidualIndicators3D result;
    result.flux_jumps = evaluate_interior_flux_jumps_3d(mesh, cell_fields);
    const auto boundary_squared = integrate_natural_residuals_3d(
        mesh, cell_fields, natural_boundary_faces, robin_boundary_faces,
        nodal_values, result.boundary_residual_faces);
    result.cell_residual_errors.resize(mesh.elements.size(), 0.0);
    result.boundary_residual_errors.resize(mesh.elements.size(), 0.0);
    result.element_errors.resize(mesh.elements.size(), 0.0);

    long double total_squared = 0.0L;
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const long double face_squared =
            static_cast<long double>(result.flux_jumps.cell_indicators[element]) *
            result.flux_jumps.cell_indicators[element];
        const long double combined_squared =
            cell_residual_squared[element] + face_squared + boundary_squared[element];
        result.cell_residual_errors[element] = checked_norm(
            cell_residual_squared[element],
            "3D cell residual indicator exceeds the finite range.");
        result.boundary_residual_errors[element] = checked_norm(
            boundary_squared[element],
            "3D boundary residual indicator exceeds the finite range.");
        result.element_errors[element] = checked_norm(
            combined_squared,
            "3D element residual indicator exceeds the finite range.");
        total_squared += combined_squared;
    }
    result.total_error = checked_norm(
        total_squared, "3D total residual indicator exceeds the finite range.");
    return result;
}


// Existing five-argument API retains its Neumann-only semantics.
DiffusionResidualIndicators3D estimate_diffusion_residual_error_3d(
    const TetrahedralMesh3D& mesh,
    const std::span<const double> nodal_values,
    const std::span<const double> element_conductivity,
    const std::function<double(double, double, double)>& source,
    const std::span<const NeumannResidualFace3D> natural_boundary_faces
) {
    return estimate_diffusion_residual_error_3d(
        mesh, nodal_values, element_conductivity, source,
        natural_boundary_faces, std::span<const RobinResidualFace3D>{});
}

// Existing four-argument API retains its cell-and-interior-only meaning.
DiffusionResidualIndicators3D estimate_diffusion_residual_error_3d(
    const TetrahedralMesh3D& mesh,
    const std::span<const double> nodal_values,
    const std::span<const double> element_conductivity,
    const std::function<double(double, double, double)>& source
) {
    return estimate_diffusion_residual_error_3d(
        mesh, nodal_values, element_conductivity, source,
        std::span<const NeumannResidualFace3D>{});
}

} // namespace nova::fem
