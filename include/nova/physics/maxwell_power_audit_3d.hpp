#pragma once
#include "nova/physics/maxwell_3d.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nova::physics {

// All powers are in watts for peak phasors and F=i*omega*mu0*J.
// Valid for homogeneous PEC and/or homogeneous natural boundaries only.
struct MaxwellPowerAudit3D {
    double source_power{};
    double operator_absorbed_power{};
    double cell_absorbed_power{};
    double residual_work{};
    double relative_power_imbalance{};
    double relative_operator_cell_mismatch{};
    double algebraic_residual_norm{};
};

// Independently compare weak source work, assembled operator loss and
// physical cell-integrated absorption. No assumption of convergence.
[[nodiscard]] inline MaxwellPowerAudit3D audit_maxwell_power_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials,
    double vacuum_wavenumber, const MaxwellSystem3D& system,
    std::span<const std::complex<double>> edge_load,
    std::span<const std::complex<double>> edge_field) {
    const auto n=system.free_edges.size(), m=system.edges.size();
    constexpr double c=299792458.0, eps0=8.8541878128e-12;
    constexpr double mu0=1.0/(eps0*c*c);
    const double omega=vacuum_wavenumber*c;
    const double denominator=2.0*omega*mu0;
    if (!(denominator>0.0) || !std::isfinite(denominator) ||
        n==0 || n>std::numeric_limits<std::size_t>::max()/2 ||
        system.real_block.rows()!=2*n || system.real_block.columns()!=2*n ||
        edge_load.size()!=m || edge_field.size()!=m)
        throw std::invalid_argument("Invalid Maxwell power audit dimensions or frequency");
    auto finite=[](std::complex<double> z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    };
    for (std::size_t i=0;i<m;++i)
        if (!finite(edge_load[i]) || !finite(edge_field[i]))
            throw std::invalid_argument("Nonfinite Maxwell power audit input");
    std::vector<double> x(2*n), b(2*n);
    std::vector<bool> free(m,false);
    for (std::size_t i=0;i<n;++i) {
        const auto j=system.free_edges[i];
        if (j>=m || (i && j<=system.free_edges[i-1]))
            throw std::invalid_argument("Invalid Maxwell free-edge ordering");
        free[j]=true;
        x[i]=edge_field[j].real(); x[i+n]=edge_field[j].imag();
        b[i]=edge_load[j].real(); b[i+n]=edge_load[j].imag();
    }
    for (std::size_t i=0;i<m;++i)
        if (!free[i] && edge_field[i]!=std::complex<double>{})
            throw std::invalid_argument("Nonzero constrained Maxwell edge field");
    const auto ax=system.real_block.multiply(x);
    MaxwellPowerAudit3D out;
    for (std::size_t i=0;i<n;++i) {
        // -Im(conj(E)*F) / (2 omega mu0) is delivered source power.
        out.source_power+=(x[i+n]*b[i]-x[i]*b[i+n])/denominator;
        out.operator_absorbed_power+=(x[i+n]*ax[i]-x[i]*ax[i+n])/denominator;
        const double rr=b[i]-ax[i], ri=b[i+n]-ax[i+n];
        out.residual_work+=(x[i+n]*rr-x[i]*ri)/denominator;
        out.algebraic_residual_norm=std::hypot(out.algebraic_residual_norm,
                                               std::hypot(rr,ri));
    }
    for (double p : integrate_maxwell_absorption_3d(
             mesh,materials,vacuum_wavenumber,system,edge_field))
        out.cell_absorbed_power+=p;
    const auto relative=[](double a,double b) {
        const double s=std::max(std::abs(a),std::abs(b));
        return s==0.0 ? 0.0 : std::abs(a/s-b/s);
    };
    out.relative_power_imbalance=relative(out.source_power,out.cell_absorbed_power);
    out.relative_operator_cell_mismatch=relative(out.operator_absorbed_power,
                                                  out.cell_absorbed_power);
    for (double v : {out.source_power,out.operator_absorbed_power,
                     out.cell_absorbed_power,out.residual_work,
                     out.relative_power_imbalance,out.relative_operator_cell_mismatch,
                     out.algebraic_residual_norm})
        if (!std::isfinite(v)) throw std::overflow_error("Maxwell power audit overflow");
    return out;
}

} // namespace nova::physics
