#include "nova/fem/nedelec_tetrahedron_3d.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nova::fem {
namespace {
double dot(const EdgeVector3D& a, const EdgeVector3D& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
EdgeVector3D cross(const EdgeVector3D& a, const EdgeVector3D& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]};
}
}

NedelecTetrahedron3D make_nedelec_tetrahedron_3d(
    const TetrahedralMesh3D& mesh, std::size_t element) {
    const auto p1 = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, element);
    const auto& nodes = mesh.elements.at(element);
    NedelecTetrahedron3D out;
    out.gradients = p1.gradients;
    out.volume = std::abs(p1.signed_jacobian_determinant)/6.0;
    if (!(out.volume > 0.0) || !std::isfinite(out.volume))
        throw std::invalid_argument("Invalid Nedelec tetrahedron volume.");
    for (const auto& g : out.gradients)
        for (const double x : g)
            if (!std::isfinite(x)) throw std::invalid_argument("Invalid Nedelec gradient.");
    for (std::size_t e=0; e<6; ++e) {
        const auto [i,j] = tetrahedron_local_edges[e];
        out.edges[e] = {std::min(nodes[i],nodes[j]), std::max(nodes[i],nodes[j])};
        out.signs[e] = nodes[i]<nodes[j] ? 1.0 : -1.0;
        out.curls[e] = cross(out.gradients[i],out.gradients[j]);
        for (double& x : out.curls[e]) {
            x *= 2.0*out.signs[e];
            if (!std::isfinite(x)) throw std::overflow_error("Nedelec curl overflow.");
        }
    }
    return out;
}

std::array<EdgeVector3D, 6> evaluate_nedelec_basis_3d(
    const NedelecTetrahedron3D& element, double r, double s, double t) {
    const auto l = evaluate_tetrahedron_p1_basis_3d(r,s,t).values;
    std::array<EdgeVector3D, 6> out{};
    for (std::size_t e=0; e<6; ++e) {
        const auto [i,j] = tetrahedron_local_edges[e];
        for (std::size_t d=0; d<3; ++d)
            out[e][d] = element.signs[e]*(l[i]*element.gradients[j][d]
                                        -l[j]*element.gradients[i][d]);
    }
    return out;
}

NedelecElementMatrices3D nedelec_element_matrices_3d(
    const NedelecTetrahedron3D& element) {
    NedelecElementMatrices3D out;
    const auto integral = [&](std::size_t a,std::size_t b) {
        return element.volume*(a==b ? 0.1 : 0.05);
    };
    const auto& g=element.gradients;
    for (std::size_t e=0; e<6; ++e) {
        const auto [i,j]=tetrahedron_local_edges[e];
        for (std::size_t f=e; f<6; ++f) {
            const auto [k,l]=tetrahedron_local_edges[f];
            const double m=element.signs[e]*element.signs[f]*(
                integral(i,k)*dot(g[j],g[l])-integral(i,l)*dot(g[j],g[k])
               -integral(j,k)*dot(g[i],g[l])+integral(j,l)*dot(g[i],g[k]));
            const double c=element.volume*dot(element.curls[e],element.curls[f]);
            if (!std::isfinite(m) || !std::isfinite(c))
                throw std::overflow_error("Nedelec element matrix overflow.");
            out.mass[e][f]=out.mass[f][e]=m;
            out.curl_curl[e][f]=out.curl_curl[f][e]=c;
        }
    }
    return out;
}

} // namespace nova::fem
