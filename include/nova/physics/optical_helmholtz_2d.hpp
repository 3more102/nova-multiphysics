#pragma once
#include "nova/fem/mesh2d.hpp"
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace nova::optics {
enum class ScalarPolarization2D { TE_Ez, TM_Hz };
struct HelmholtzMaterial2D {
    std::complex<double> epsilon_r{1.0};
    std::complex<double> mu_r{1.0};
};
struct HelmholtzDirichlet2D {
    std::size_t node{};
    std::complex<double> value{};
};
// a du/dn + alpha*u = load; unassigned exterior edges have zero flux.
struct HelmholtzRobin2D {
    fem::BoundaryEdge2D edge{};
    std::complex<double> alpha{};
    std::complex<double> load{};
};
struct Helmholtz2DResult {
    std::vector<std::complex<double>> field;
    std::vector<std::array<std::complex<double>, 2>> cell_gradient;
    // Nonnegative normalized loss integrals, not watts.
    std::vector<double> cell_loss;
    double absolute_residual{};
    double relative_backward_error{};
    // Signed normalized weak-form power terms (not watts).
    double excitation_input{}; // volume and Robin loads
    double dirichlet_input{};  // reactions on fixed nodes
    double robin_outflow{};    // homogeneous Robin impedance
    double material_loss{};
    double power_balance_residual{}; // inputs - outflow - loss
};
// Experimental bounded dense P1 reference. SI mesh coordinates and k0 in rad/m.
// TE: u=Ez, a=1/mu_r, b=epsilon_r; TM: u=Hz, a=1/epsilon_r, b=mu_r.
// Solves -div(a grad u)-k0^2 b u=f, exp(-i omega t).
// One material and optionally one constant source per triangle.
[[nodiscard]] Helmholtz2DResult solve_helmholtz_2d(
    const fem::TriangleMesh2D& mesh, double k0, ScalarPolarization2D polarization,
    std::span<const HelmholtzMaterial2D> materials,
    std::span<const HelmholtzDirichlet2D> prescribed = {},
    std::span<const HelmholtzRobin2D> robin = {},
    std::span<const std::complex<double>> sources = {});
} // namespace nova::optics
