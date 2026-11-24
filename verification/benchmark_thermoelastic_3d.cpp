// Matched workload against the pre-coupling SDK. Build the same file with
// NOVA_BENCH_REFERENCE and the baseline headers/library to time the old path.
#include "nova/fem/p1_flux_3d.hpp"
#include "nova/fem/transient_heat_3d.hpp"
#include "nova/io/gmsh_tetrahedral.hpp"
#include "nova/io/vtk_legacy_3d.hpp"
#include "nova/physics/linear_elasticity_3d.hpp"
#ifndef NOVA_BENCH_REFERENCE
#include "nova/physics/thermoelasticity_3d.hpp"
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

int main(int argc, char** argv) {
    using namespace nova;
    using Clock = std::chrono::steady_clock;
    const auto seconds = [](auto a, auto b) { return std::chrono::duration<double>(b-a).count(); };
    try {
        if (argc != 3) throw std::invalid_argument("Usage: benchmark mesh.msh final_heat.vtk");
        const auto start = Clock::now();
        std::ifstream input(argv[1]);
        auto imported = io::read_gmsh_tetrahedra(input);
        const auto mesh_done = Clock::now();
        fem::TransientHeat3DProblem heat;
        heat.diffusion.mesh = std::move(imported.mesh);
        const auto& mesh = heat.diffusion.mesh;
        heat.diffusion.cell_diffusivity.assign(mesh.elements.size(), 2.);
        heat.diffusion.source = [](double,double,double) { return -4.; };
        heat.cell_volumetric_capacity.assign(mesh.elements.size(),2.);
        heat.initial_temperature.assign(mesh.nodes.size(),310.);
        std::vector<physics::DisplacementConstraint3D> fixed;
        for (std::size_t n = 0; n < mesh.nodes.size(); ++n) {
            if (mesh.nodes[n].x == 0) fixed.push_back({n,{}});
        }
        std::vector<physics::IsotropicElasticMaterial3D> material;
        for (const auto tag : imported.volume_physical_tags) material.push_back({tag == 11 ? 1000. : 2000.,0.2});
#ifdef NOVA_BENCH_REFERENCE
        physics::LinearElasticity3DProblem legacy;
        legacy.mesh = mesh;
        legacy.cell_materials = std::move(material);
        legacy.prescribed_displacements = std::move(fixed);
        legacy.cell_thermal_expansion_per_kelvin.assign(mesh.elements.size(),.001);
        legacy.cell_temperature_change_kelvin.resize(mesh.elements.size());
#else
        physics::PreparedElasticity3D prepared(mesh,material,fixed);
        const std::vector<physics::ThermalExpansion3D> expansion(mesh.elements.size(),{.001,300.});
#endif
        const auto setup_done = Clock::now();
        double responses = 0, energy_sum = 0, displacement_norm = 0, stress_norm = 0, temperature_error = 0;
        const linalg::SolverOptions linear{2000,1e-10,1e-14,linalg::PreconditionerType::Jacobi};
        fem::TransientHeat3DOptions options;
        options.scheme = fem::TransientHeat3DScheme::CrankNicolson;
        options.time_step = .125;
        options.steps = 5;
        options.retain_step_solutions = false;
        options.solver_options = linear;
        options.on_step = [&](std::size_t, double time, const std::vector<double>& t) {
            const auto before = Clock::now();
#ifdef NOVA_BENCH_REFERENCE
            for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
                double change = 0;
                for (const auto n : mesh.elements[e]) change += .25*(t[n]-300.);
                legacy.cell_temperature_change_kelvin[e] = change;
            }
            const auto result = physics::solve_linear_elasticity_3d(legacy,linear);
#else
            const auto result = prepared.solve_temperature(t,expansion,linear);
#endif
            if (!result.linear_solver.converged) throw std::runtime_error(result.linear_solver.message);
            energy_sum += result.strain_energy_joules;
            for (const auto& u : result.displacement) for (const double v : u) displacement_norm += v*v;
            for (const auto& s : result.cell_stress) for (const double v : s) stress_norm += v*v;
            for (const double v : t) temperature_error = std::max(temperature_error,std::abs(v-(310.-2.*time)));
            responses += seconds(before,Clock::now());
        };
        const auto result = fem::solve_transient_heat_3d(heat,options);
        const auto solved = Clock::now();
        auto flux = fem::evaluate_p1_tetrahedral_flux_3d(heat.diffusion,result.temperature);
        std::ofstream output(argv[2]);
        io::write_legacy_vtk_conduction_3d(output,mesh,result.temperature,imported.volume_physical_tags,flux);
        output.close();
        if (!output || temperature_error > 5e-8) throw std::runtime_error("Benchmark failed field/output check.");
        const auto finished = Clock::now();
        std::cout << std::setprecision(15)
                  << "mesh_s " << seconds(start,mesh_done) << " setup_s " << seconds(mesh_done,setup_done)
                  << " heat_s " << seconds(setup_done,solved)-responses << " response_s " << responses
                  << " export_s " << seconds(solved,finished) << " total_s " << seconds(start,finished)
                  << " energy " << energy_sum << " displacement_norm " << std::sqrt(displacement_norm)
                  << " stress_norm " << std::sqrt(stress_norm) << " temperature_error " << temperature_error;
#ifndef NOVA_BENCH_REFERENCE
        const auto profile = prepared.timings();
        std::cout << " topology_s " << profile.topology << " element_setup_s " << profile.element_setup
                  << " stiffness_s " << profile.stiffness_assembly << " thermal_rhs_s " << profile.thermal_load
                  << " cg_s " << profile.linear_solve << " recovery_s " << profile.recovery
                  << " elastic_nnz " << prepared.stiffness_nonzeros();
#endif
        std::cout << '\n';
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
