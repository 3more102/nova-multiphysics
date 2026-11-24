#pragma once

#include "nova/physics/maxwell_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <set>
#include <stdexcept>
#include <vector>

namespace nova::physics {

// Physical SI coordinates, cell index; returns the vector weak source F,
// not the current J. For exp(-i omega t), F = i omega mu0 J.
using MaxwellSpatialSource3D = std::function<ComplexVector3D(
    std::size_t, const std::array<double,3>&)>;

// Four-point Gauss-Legendre in each Duffy coordinate (64 samples/cell).
// Exact through degree five for polynomial integrands on affine tetrahedra,
// hence through degree four for F times the lowest-order Nedelec basis.
// Returns full globally oriented edge loads, including constrained edges.
[[nodiscard]] inline std::vector<std::complex<double>>
assemble_maxwell_spatial_load_3d(
    const fem::TetrahedralMesh3D& mesh,
    const MaxwellSpatialSource3D& source,
    const MaxwellSystem3D& system) {
    mesh.validate();
    if (!source) throw std::invalid_argument("Missing Maxwell spatial source.");
    std::set<fem::Edge3D> expected;
    for (const auto& cell : mesh.elements)
        for (const auto& [i,j] : fem::tetrahedron_local_edges)
            expected.insert({std::min(cell[i],cell[j]),std::max(cell[i],cell[j])});
    if (expected.size()!=system.edges.size() ||
        !std::equal(expected.begin(),expected.end(),system.edges.begin()))
        throw std::invalid_argument("Maxwell spatial source edge topology mismatch.");

    constexpr std::array<double,4> q{
        0.06943184420297371,0.33000947820757187,
        0.66999052179242813,0.93056815579702629};
    constexpr std::array<double,4> w{
        0.17392742256872693,0.32607257743127307,
        0.32607257743127307,0.17392742256872693};
    const auto finite=[](std::complex<double> z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    };
    std::vector<std::complex<double>> load(system.edges.size());
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        const auto cell=fem::make_nedelec_tetrahedron_3d(mesh,c);
        const auto& ids=mesh.elements[c];
        std::array<std::size_t,6> global{};
        for (std::size_t e=0;e<6;++e)
            global[e]=static_cast<std::size_t>(std::lower_bound(
                system.edges.begin(),system.edges.end(),cell.edges[e])-system.edges.begin());
        for (std::size_t i=0;i<4;++i)
            for (std::size_t j=0;j<4;++j)
                for (std::size_t k=0;k<4;++k) {
                    const double u=q[i], v=q[j], h=q[k];
                    const double r=u, s=(1-u)*v, t=(1-u)*(1-v)*h;
                    const std::array<double,4> l{1-r-s-t,r,s,t};
                    std::array<double,3> p{};
                    for (std::size_t n=0;n<4;++n) {
                        const auto& x=mesh.nodes[ids[n]];
                        p[0]+=l[n]*x.x; p[1]+=l[n]*x.y; p[2]+=l[n]*x.z;
                    }
                    for (double x:p)
                        if (!std::isfinite(x))
                            throw std::overflow_error("Maxwell quadrature coordinate overflow.");
                    const auto f=source(c,p);
                    for (auto z:f)
                        if (!finite(z))
                            throw std::invalid_argument("Nonfinite Maxwell spatial source.");
                    const auto basis=fem::evaluate_nedelec_basis_3d(cell,r,s,t);
                    const double weight=cell.volume*(6*w[i]*w[j]*w[k]
                        *(1-u)*(1-u)*(1-v));
                    for (std::size_t e=0;e<6;++e) {
                        auto& value=load[global[e]];
                        for (std::size_t d=0;d<3;++d)
                            value+=(weight*basis[e][d])*f[d];
                        if (!finite(value))
                            throw std::overflow_error("Maxwell spatial load overflow.");
                    }
                }
    }
    return load;
}
} // namespace nova::physics
