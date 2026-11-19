#include "nova/fem/adaptive_heat_3d.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace {
using namespace nova::fem;

constexpr auto usage = "Usage: nova_anisotropic_heat_3d [be|cn [subdivisions [atol [--csv]]]]";

template<class Number>
Number parse_number(std::string_view input) {
    Number value{};
    const auto end = input.data() + input.size();
    const auto parsed = std::from_chars(input.data(), end, value);
    if (parsed.ec != std::errc{} || parsed.ptr != end) throw std::invalid_argument(usage);
    return value;
}

TransientHeat3DProblem cooling_cube(std::size_t n) {
    auto box = make_box_mesh_3d(0,1,0,1,0,1, n,n,n);
    TransientHeat3DProblem problem;
    problem.mesh = std::move(box.mesh);
    problem.volumetric_capacity = [](double,double,double) { return 1.0; };
    problem.conductivity_tensor = [](double x,double,double) {
        const double s = 1 + 0.5*x;
        return SymmetricTensor3D{2*s, 1.5*s, s, 0.4*s, -0.25*s, 0.3*s};
    };
    problem.initial_temperature = [](double x,double y,double z) {
        return 64*x*(1-x)*y*(1-y)*z*(1-z)*(1+0.5*x);
    };
    problem.source = [](double,double,double z,double t) { return (1+t)*(1+z); };
    for (const auto& face : box.boundary_faces) {
        problem.dirichlet_faces.push_back({face, [](double,double,double,double) { return 0.0; }});
    }
    return problem;
}
} // namespace

int main(int argc, char** argv) {
    try {
        const std::string_view scheme = argc > 1 ? argv[1] : "cn";
        const auto n = argc > 2 ? parse_number<std::size_t>(argv[2]) : 4U;
        const double tolerance = argc > 3 ? parse_number<double>(argv[3]) : 1e-5;
        const bool csv = argc == 5 && std::string_view(argv[4]) == "--csv";
        if ((argc > 4 && !csv) || (scheme != "be" && scheme != "cn") ||
            n < 2 || !(tolerance > 0) || !std::isfinite(tolerance)) {
            throw std::invalid_argument(usage);
        }
        const auto problem = cooling_cube(n);
        AdaptiveHeat3DOptions options;
        options.scheme = scheme == "be" ? TransientTimeScheme::BackwardEuler : TransientTimeScheme::CrankNicolson;
        options.end_time = options.initial_step = options.max_step = 0.02;
        options.min_step = 1e-10;
        options.absolute_tolerance = tolerance;
        options.relative_tolerance = 0;
        options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
        options.solver_options.relative_tolerance = 1e-13;
        options.solver_options.absolute_tolerance = 1e-15;
        const auto result = solve_adaptive_heat_3d(problem, options);
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
        if (csv) {
            for (const auto& trial : result.attempts) {
                std::cout << "#trial," << trial.start_time << ',' << trial.attempted_step << ','
                          << trial.normalized_error << ',' << trial.accepted << '\n';
            }
            std::cout << "x,y,z,temperature\n";
            for (std::size_t i = 0; i < result.temperature.size(); ++i) {
                const auto& x = problem.mesh.nodes[i];
                std::cout << x.x << ',' << x.y << ',' << x.z << ',' << result.temperature[i] << '\n';
            }
        } else {
            double max_imbalance = 0;
            for (const auto& d : result.diagnostics) {
                max_imbalance = std::max(max_imbalance, std::abs(d.balance_residual));
            }
            std::cout << "tetrahedra=" << problem.mesh.elements.size() << " dofs=" << result.temperature.size()
                      << " final_time=" << result.final_time << " accepted=" << result.accepted_steps
                      << " rejected=" << result.rejected_steps << " max_energy_imbalance=" << max_imbalance << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
