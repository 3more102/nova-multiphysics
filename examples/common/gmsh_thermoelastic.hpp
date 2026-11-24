#pragma once

#include "nova/io/gmsh_tetrahedral.hpp"
#include "nova/io/vtk_legacy_3d.hpp"
#include "nova/physics/thermoelasticity_3d.hpp"

#include <charconv>
#include <cmath>
#include <iomanip>
#include <map>
#include <optional>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace nova::examples {

template<class T>
T number(std::string_view token) {
    T value{};
    const auto end = token.data() + token.size();
    const auto parsed = std::from_chars(token.data(), end, value);
    if (parsed.ec != std::errc{} || parsed.ptr != end || !std::isfinite(value)) {
        throw std::invalid_argument("Invalid finite number: " + std::string(token));
    }
    return value;
}

using io::physical_tag;

struct ElasticRegion {
    physics::IsotropicElasticMaterial3D material;
    physics::ThermalExpansion3D expansion;
};

struct TaggedMechanics {
    std::map<int, ElasticRegion> volumes;
    std::set<int> clamps;
    std::map<std::pair<int, std::size_t>, double> components;
    std::optional<linalg::SolverOptions> solver;

    bool enabled() const { return !volumes.empty() || !clamps.empty() || !components.empty() || solver.has_value(); }

    template<class Next>
    bool read(std::string_view option, Next& next, const io::GmshTetrahedralMesh3D& mesh) {
        if (option == "--elastic") {
            const auto tag = physical_tag(next(), 3, mesh);
            const ElasticRegion law{{number<double>(next()), number<double>(next())},
                                    {number<double>(next()), number<double>(next())}};
            if (!volumes.emplace(tag, law).second) throw std::invalid_argument("Duplicate elastic volume.");
        } else if (option == "--clamp") {
            if (!clamps.insert(physical_tag(next(), 2, mesh)).second) {
                throw std::invalid_argument("Duplicate clamped surface group.");
            }
        } else if (option == "--fix-component") {
            const auto tag = physical_tag(next(), 2, mesh);
            const auto axis = next();
            const auto index = std::string_view("xyz").find(axis);
            if (axis.size() != 1 || index == std::string_view::npos) {
                throw std::invalid_argument("Displacement component must be x, y or z.");
            }
            if (!components.emplace(std::make_pair(tag, index), number<double>(next())).second) {
                throw std::invalid_argument("Duplicate prescribed surface component.");
            }
        } else if (option == "--elastic-solver") {
            if (solver) throw std::invalid_argument("Duplicate elastic solver settings.");
            const auto kind = next();
            linalg::SolverOptions options;
            if (kind == "jacobi") options.preconditioner = linalg::PreconditionerType::Jacobi;
            else if (kind == "ic0") options.preconditioner = linalg::PreconditionerType::IncompleteCholesky0;
            else throw std::invalid_argument("Elastic preconditioner must be jacobi or ic0.");
            options.relative_tolerance = number<double>(next());
            options.max_iterations = number<std::size_t>(next());
            if (!(options.relative_tolerance > 0 && options.relative_tolerance < 1) || options.max_iterations == 0) {
                throw std::invalid_argument("Elastic solver requires 0 < rtol < 1 and a positive iteration limit.");
            }
            solver = options;
        } else {
            return false;
        }
        return true;
    }
};

struct TaggedThermoelasticStudy {
    physics::PreparedElasticity3D stiffness;
    std::vector<physics::ThermalExpansion3D> expansion;
    linalg::SolverOptions options;

    physics::LinearElasticity3DResult solve(const std::vector<double>& temperature) {
        auto result = stiffness.solve_temperature(temperature, expansion, options);
        if (!result.linear_solver.converged) {
            throw std::runtime_error("Thermoelastic CG failed: " + result.linear_solver.message);
        }
        return result;
    }
};

inline std::optional<TaggedThermoelasticStudy> prepare_mechanics(
    const io::GmshTetrahedralMesh3D& mesh, const TaggedMechanics& settings
) {
    if (!settings.enabled()) return std::nullopt;
    const std::set<int> volume_tags(mesh.volume_physical_tags.begin(), mesh.volume_physical_tags.end());
    const std::set<int> face_tags(mesh.exterior_physical_tags.begin(), mesh.exterior_physical_tags.end());
    for (const auto& [tag, law] : settings.volumes) {
        if (!volume_tags.contains(tag)) throw std::invalid_argument("Unknown elastic volume tag " + std::to_string(tag));
    }
    for (const auto tag : settings.clamps) {
        if (!face_tags.contains(tag)) throw std::invalid_argument("Unknown clamp surface tag " + std::to_string(tag));
    }
    for (const auto& entry : settings.components) {
        if (!face_tags.contains(entry.first.first)) throw std::invalid_argument("Unknown displacement surface tag.");
    }
    std::vector<physics::IsotropicElasticMaterial3D> materials;
    std::vector<physics::ThermalExpansion3D> expansion;
    materials.reserve(mesh.volume_physical_tags.size());
    expansion.reserve(mesh.volume_physical_tags.size());
    for (const auto tag : mesh.volume_physical_tags) {
        const auto region = settings.volumes.find(tag);
        if (region == settings.volumes.end()) throw std::invalid_argument("Missing elastic material for volume " + std::to_string(tag));
        materials.push_back(region->second.material);
        expansion.push_back(region->second.expansion);
    }
    std::map<std::pair<std::size_t, std::size_t>, double> nodal;
    for (std::size_t f = 0; f < mesh.exterior_faces.size(); ++f) {
        const auto tag = mesh.exterior_physical_tags[f];
        const bool clamped = settings.clamps.contains(tag);
        for (std::size_t i = 0; i < 3; ++i) {
            const auto bc = settings.components.find({tag,i});
            if (!clamped && bc == settings.components.end()) continue;
            const auto value = bc == settings.components.end() ? 0.0 : bc->second;
            for (const auto node : mesh.exterior_faces[f]) {
                const auto [it, added] = nodal.emplace(std::make_pair(node,i), value);
                if ((clamped && value != 0.0) || (!added && it->second != value)) {
                    throw std::invalid_argument("Conflicting prescribed displacement on a shared vertex.");
                }
            }
        }
    }
    std::vector<physics::DisplacementComponent3D> constraints;
    constraints.reserve(nodal.size());
    for (const auto& [dof,value] : nodal) constraints.push_back({dof.first, dof.second, value});
    linalg::SolverOptions defaults;
    defaults.relative_tolerance = 1e-11;
    defaults.preconditioner = linalg::PreconditionerType::Jacobi;
    return TaggedThermoelasticStudy{physics::PreparedElasticity3D(mesh.mesh, materials, {}, {}, {}, constraints),
                                    std::move(expansion), settings.solver.value_or(defaults)};
}

inline void write_heat_fields(
    std::ostream& out, const fem::TetrahedralMesh3D& mesh, const std::vector<double>& temperature,
    const std::vector<int>& tags, const std::vector<fem::P1TetrahedralFlux3D>& flux,
    const std::vector<double>& generation, const std::optional<physics::LinearElasticity3DResult>& response
) {
    if (response) io::write_legacy_vtk_thermoelastic_3d(out, mesh, temperature, tags, flux, *response, generation);
    else io::write_legacy_vtk_conduction_3d(out, mesh, temperature, tags, flux, generation);
}

inline void report_mechanics(std::ostream& out, const physics::LinearElasticity3DResult& result) {
    out << std::setprecision(17) << "elastic_iterations " << result.linear_solver.iterations
        << " elastic_residual " << result.linear_solver.residual_norm
        << " elastic_strain_energy_J " << result.strain_energy_joules << '\n';
}

} // namespace nova::examples
