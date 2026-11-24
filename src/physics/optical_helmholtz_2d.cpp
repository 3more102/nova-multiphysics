#include "nova/physics/optical_helmholtz_2d.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace nova::optics {
namespace {
using C = std::complex<double>;
using Edge = std::pair<std::size_t, std::size_t>;
bool finite(C z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }
Edge edge_key(std::size_t i, std::size_t j) { return std::minmax(i,j); }
struct Cell {
    double area;
    std::array<std::array<double,2>,3> grad;
    C a, b;
};
Cell geometry(const fem::TriangleMesh2D& mesh, const fem::Triangle2D& e,
              C a, C b) {
    const auto p = mesh.nodes[e[0]], q = mesh.nodes[e[1]], r = mesh.nodes[e[2]];
    const double x = q.x-p.x, y = q.y-p.y, s = r.x-p.x, t = r.y-p.y;
    const double d = x*t-y*s;
    return {std::abs(d)/2, {{{(y-t)/d,(s-x)/d},{t/d,-s/d},{-y/d,x/d}}}, a,b};
}
std::vector<C> dense_solve(std::vector<C> a, std::vector<C> b) {
    const auto n = b.size();
    if (n == 0) return {};
    // Equilibrate rows before partial pivoting; reject unresolved pivots.
    for (std::size_t i=0;i<n;++i) {
        double scale=0;
        for (std::size_t j=0;j<n;++j) scale=std::max(scale,std::abs(a[i*n+j]));
        if (!(scale>0) || !std::isfinite(scale)) throw std::runtime_error("Invalid Helmholtz row");
        for (std::size_t j=0;j<n;++j) a[i*n+j]/=scale;
        b[i]/=scale;
        if (!finite(b[i])) throw std::runtime_error("Nonfinite scaled Helmholtz load");
    }
    constexpr double threshold=64*std::numeric_limits<double>::epsilon();
    for (std::size_t k=0;k<n;++k) {
        std::size_t p=k;
        for (std::size_t i=k+1;i<n;++i)
            if (std::abs(a[i*n+k])>std::abs(a[p*n+k])) p=i;
        if (!(std::abs(a[p*n+k])>threshold) || !finite(a[p*n+k]))
            throw std::runtime_error("Singular or numerically unresolved Helmholtz pivot");
        if (p!=k) {
            for (std::size_t j=k;j<n;++j) std::swap(a[p*n+j],a[k*n+j]);
            std::swap(b[p],b[k]);
        }
        for (std::size_t i=k+1;i<n;++i) {
            const C f=a[i*n+k]/a[k*n+k];
            a[i*n+k]=0;
            if (f==C{}) continue;
            for (std::size_t j=k+1;j<n;++j) a[i*n+j]-=f*a[k*n+j];
            b[i]-=f*b[k];
        }
    }
    for (std::size_t i=n;i-- > 0;) {
        for (std::size_t j=i+1;j<n;++j) b[i]-=a[i*n+j]*b[j];
        b[i]/=a[i*n+i];
        if (!finite(b[i])) throw std::runtime_error("Nonfinite Helmholtz solution");
    }
    return b;
}
} // namespace

Helmholtz2DResult solve_helmholtz_2d(
    const fem::TriangleMesh2D& mesh, double k0, ScalarPolarization2D polarization,
    std::span<const HelmholtzMaterial2D> materials,
    std::span<const HelmholtzDirichlet2D> prescribed,
    std::span<const HelmholtzRobin2D> robin, std::span<const C> sources) {
    const auto n=mesh.nodes.size(), m=mesh.elements.size();
    if (n>1024) throw std::length_error("Dense Helmholtz reference limited to 1024 nodes");
    mesh.validate();
    const double k2=k0*k0;
    if (!(k0>0) || !(k2>0) || !std::isfinite(k2) || materials.size()!=m ||
        (!sources.empty() && sources.size()!=m))
        throw std::invalid_argument("Invalid Helmholtz frequency or cell data");
    const bool te=polarization==ScalarPolarization2D::TE_Ez;
    if (!te && polarization!=ScalarPolarization2D::TM_Hz)
        throw std::invalid_argument("Invalid scalar polarization");
    std::vector<C> A(n*n), rhs(n), field(n);
    std::vector<bool> fixed(n,false);
    std::vector<Cell> cells;
    std::map<Edge,unsigned> incidence;
    for (std::size_t c=0;c<m;++c) {
        const auto mat=materials[c];
        if (!finite(mat.epsilon_r) || !finite(mat.mu_r) ||
            mat.epsilon_r.imag()<0 || mat.mu_r.imag()<0)
            throw std::invalid_argument("Expected passive finite scalar material");
        const C denominator=te?mat.mu_r:mat.epsilon_r;
        if (denominator==C{}) throw std::invalid_argument("Zero inverse-weight material");
        const C a=1.0/denominator, b=te?mat.epsilon_r:mat.mu_r;
        if (!finite(a) || a==C{}) throw std::invalid_argument("Unrepresentable inverse material");
        const auto& e=mesh.elements[c];
        const Cell cell=geometry(mesh,e,a,b);
        cells.push_back(cell);
        const C f=sources.empty()?C{}:sources[c];
        if (!finite(f)) throw std::invalid_argument("Nonfinite Helmholtz source");
        for (std::size_t i=0;i<3;++i) {
            ++incidence[edge_key(e[i],e[(i+1)%3])];
            rhs[e[i]]+=cell.area*f/3.0;
            for (std::size_t j=0;j<3;++j) {
                const double g=cell.grad[i][0]*cell.grad[j][0]+cell.grad[i][1]*cell.grad[j][1];
                A[e[i]*n+e[j]]+=a*(cell.area*g)-k2*b*(cell.area*(i==j?2.0:1.0)/12.0);
            }
        }
    }
    std::set<Edge> assigned;
    std::vector<bool> exterior(n,false);
    for (const auto& [edge,count]:incidence) if (count==1) {
        exterior[edge.first]=true; exterior[edge.second]=true;
    }
    for (const auto& bc:prescribed) {
        if (bc.node>=n || !exterior[bc.node] || fixed[bc.node] || !finite(bc.value))
            throw std::invalid_argument("Invalid or repeated exterior Dirichlet node");
        fixed[bc.node]=true; field[bc.node]=bc.value;
    }
    for (const auto& bc:robin) {
        const auto edge=edge_key(bc.edge[0],bc.edge[1]);
        if (!incidence.contains(edge) || incidence.at(edge)!=1 ||
            !assigned.insert(edge).second || !finite(bc.alpha) || !finite(bc.load))
            throw std::invalid_argument("Invalid or repeated exterior Robin edge");
        if (fixed[edge.first] && fixed[edge.second])
            throw std::invalid_argument("Robin edge fully constrained by Dirichlet data");
        const auto p=mesh.nodes[edge.first], q=mesh.nodes[edge.second];
        const double length=std::hypot(p.x-q.x,p.y-q.y);
        for (std::size_t i:bc.edge) {
            rhs[i]+=bc.load*length/2.0;
            for (std::size_t j:bc.edge) A[i*n+j]+=bc.alpha*length*(i==j?2.0:1.0)/6.0;
        }
    }
    for (C v:A) if (!finite(v)) throw std::runtime_error("Nonfinite Helmholtz matrix");
    for (C v:rhs) if (!finite(v)) throw std::runtime_error("Nonfinite Helmholtz load");
    std::vector<std::size_t> free;
    for (std::size_t i=0;i<n;++i) if (!fixed[i]) free.push_back(i);
    const auto nf=free.size();
    std::vector<C> reduced(nf*nf), load(nf);
    for (std::size_t i=0;i<nf;++i) {
        const auto row=free[i];
        load[i]=rhs[row];
        for (std::size_t j=0;j<n;++j) if (fixed[j]) load[i]-=A[row*n+j]*field[j];
        for (std::size_t j=0;j<nf;++j) reduced[i*nf+j]=A[row*n+free[j]];
    }
    const auto solution=dense_solve(std::move(reduced),std::move(load));
    for (std::size_t i=0;i<nf;++i) field[free[i]]=solution[i];
    Helmholtz2DResult result;
    double anorm=0, xnorm=0, bnorm=0;
    for (C u:field) xnorm=std::max(xnorm,std::abs(u));
    for (auto i:free) {
        C residual=-rhs[i];
        double rownorm=0;
        for (std::size_t j=0;j<n;++j) {
            residual+=A[i*n+j]*field[j]; rownorm+=std::abs(A[i*n+j]);
        }
        if (!finite(residual) || !std::isfinite(rownorm)) throw std::runtime_error("Nonfinite Helmholtz residual");
        result.absolute_residual=std::max(result.absolute_residual,std::abs(residual));
        anorm=std::max(anorm,rownorm); bnorm=std::max(bnorm,std::abs(rhs[i]));
    }
    const double scale=anorm*xnorm+bnorm;
    if (!std::isfinite(scale)) throw std::runtime_error("Helmholtz residual scale overflow");
    result.relative_backward_error=scale>0?result.absolute_residual/scale:result.absolute_residual;
    if (result.relative_backward_error>1e-10) throw std::runtime_error("Helmholtz residual exceeds reference tolerance");
    for (std::size_t c=0;c<m;++c) {
        const auto& e=mesh.elements[c]; const auto& cell=cells[c];
        std::array<C,2> grad{};
        double diagonal=0; C sum{};
        for (std::size_t i=0;i<3;++i) {
            const C u=field[e[i]]; diagonal+=std::norm(u); sum+=u;
            for (std::size_t d=0;d<2;++d) grad[d]+=u*cell.grad[i][d];
        }
        const double integral=cell.area*(diagonal+std::norm(sum))/12.0;
        const double loss=cell.area*(-cell.a.imag())*(std::norm(grad[0])+std::norm(grad[1]))
                        +k2*cell.b.imag()*integral;
        if (!finite(grad[0]) || !finite(grad[1]) || !std::isfinite(loss) || loss<0)
            throw std::runtime_error("Invalid reconstructed Helmholtz loss");
        result.cell_gradient.push_back(grad);
        result.cell_loss.push_back(loss);
    }
    // Imaginary part of the discrete weak form, including constrained rows.
    // For exp(-i omega t), passive material and outgoing impedance have
    // negative imaginary quadratic forms.
    for (std::size_t i=0;i<n;++i) {
        result.excitation_input-=std::imag(std::conj(field[i])*rhs[i]);
        if (fixed[i]) {
            C reaction=-rhs[i];
            for (std::size_t j=0;j<n;++j) reaction+=A[i*n+j]*field[j];
            result.dirichlet_input-=std::imag(std::conj(field[i])*reaction);
        }
    }
    for (const auto& bc:robin) {
        const auto p=mesh.nodes[bc.edge[0]], q=mesh.nodes[bc.edge[1]];
        const C u=field[bc.edge[0]], v=field[bc.edge[1]];
        const double integral=std::hypot(p.x-q.x,p.y-q.y)*
            (std::norm(u)+std::norm(v)+std::real(std::conj(u)*v))/3.0;
        result.robin_outflow-=bc.alpha.imag()*integral;
    }
    for (double loss:result.cell_loss) result.material_loss+=loss;
    result.power_balance_residual=result.excitation_input+result.dirichlet_input
        -result.robin_outflow-result.material_loss;
    if (!std::isfinite(result.excitation_input) ||
        !std::isfinite(result.dirichlet_input) ||
        !std::isfinite(result.robin_outflow) ||
        !std::isfinite(result.material_loss) ||
        !std::isfinite(result.power_balance_residual))
        throw std::runtime_error("Nonfinite Helmholtz power audit");
    result.field=std::move(field);
    return result;
}
} // namespace nova::optics
