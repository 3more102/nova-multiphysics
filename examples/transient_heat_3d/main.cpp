#include "nova/fem/transient_heat_3d.hpp"
#include "nova/io/vtk_legacy_3d.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace {
std::size_t positive_count(std::string_view arg) {
    std::size_t value{};
    const auto end = arg.data() + arg.size();
    const auto parsed = std::from_chars(arg.data(), end, value);
    if (parsed.ec != std::errc{} || parsed.ptr != end || value == 0) {
        throw std::invalid_argument("Mesh subdivisions and step count must be positive integers.");
    }
    return value;
}

double spatial_mode(double x, double y, double z) {
    constexpr double pi = std::numbers::pi;
    return std::sin(pi*x) * std::sin(pi*y) * std::sin(pi*z);
}
} // namespace

int main(int argc, char** argv) {
    using namespace nova::fem;
    try {
        const bool csv = argc == 4 && std::string_view(argv[3]) == "--csv";
        const bool vtk = argc == 5 && std::string_view(argv[3]) == "--vtk";
        if (argc > 5 || (argc > 3 && !csv && !vtk)) {
            throw std::invalid_argument("Usage: nova_transient_heat_3d [subdivisions [steps [--csv | --vtk output.vtk]]]");
        }
        const std::size_t n = argc > 1 ? positive_count(argv[1]) : 4;
        if (n < 2) throw std::invalid_argument("The cooling cube needs at least two subdivisions.");
        const auto steps = argc > 2 ? positive_count(argv[2]) : 10;
        auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, n, n, n);
        TransientHeat3DProblem problem;
        problem.mesh = std::move(box.mesh);
        problem.initial_temperature = spatial_mode;
        problem.volumetric_capacity = [](double, double, double) { return 1.0; };
        problem.conductivity = problem.volumetric_capacity;
        problem.source = [](double, double, double, double) { return 0.0; };
        for (const auto& face : box.boundary_faces) {
            problem.dirichlet_faces.push_back({face, problem.source});
        }
        TransientHeat3DOptions options;
        options.scheme = TransientTimeScheme::CrankNicolson;
        options.steps = steps;
        options.time_step = 0.02 / static_cast<double>(steps);
        options.solver_options.relative_tolerance = 1.0e-12;
        options.solver_options.absolute_tolerance = 1.0e-15;
        options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;

        const auto start = std::chrono::steady_clock::now();
        const auto result = solve_transient_heat_3d(problem, options);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
        if (argc == 4) {
            std::cout << "x,y,z,temperature\n";
            for (std::size_t i = 0; i < result.temperature.size(); ++i) {
                const auto& point = problem.mesh.nodes[i];
                std::cout << point.x << ',' << point.y << ',' << point.z << ',' << result.temperature[i] << '\n';
            }
            return 0;
        }

        if (vtk) {
            std::ofstream output(argv[4]);
            if (!output) throw std::runtime_error("Cannot open VTK output file.");
            nova::io::write_legacy_vtk_scalar_3d(output, problem.mesh, result.temperature);
            output.close();
            if (!output) throw std::runtime_error("Failed to close VTK output file.");
        }

        const double decay = std::exp(-3.0 * std::numbers::pi * std::numbers::pi * result.final_time);
        double max_error = 0.0, max_balance = 0.0;
        std::size_t iterations = 0;
        for (std::size_t i = 0; i < result.temperature.size(); ++i) {
            const auto& point = problem.mesh.nodes[i];
            const double exact = decay * spatial_mode(point.x, point.y, point.z);
            max_error = std::max(max_error, std::abs(result.temperature[i] - exact));
        }
        for (const auto& d : result.diagnostics) max_balance = std::max(max_balance, std::abs(d.balance_residual));
        for (const auto& solver : result.linear_solvers) iterations += solver.iterations;
        std::cout << "tetrahedra=" << problem.mesh.elements.size() << " dofs=" << result.temperature.size()
                  << " steps=" << steps << " cg_iterations=" << iterations << "\nmax_nodal_error=" << max_error
                  << " max_step_energy_imbalance=" << max_balance << "\nsolve_seconds=" << seconds << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
