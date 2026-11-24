#include "nova/physics/optical_helmholtz_1d.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace nova::optics {
namespace {
using C = std::complex<double>;
constexpr C I{0.0, 1.0};
bool finite(C z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }
void check(C z) {
    if (!finite(z)) throw std::overflow_error("nonfinite Helmholtz coefficient or result");
}
void positive(double x) {
    if (!std::isfinite(x) || x <= 0.0) throw std::invalid_argument("expected finite positive input");
}
// Tridiagonal elimination with adjacent partial pivoting. lower stores
// the second superdiagonal after each elimination step.
std::vector<C> solve(std::vector<C> lower, std::vector<C> diag,
                     std::vector<C> upper, std::vector<C> rhs) {
    const auto n = diag.size();
    for (std::size_t j = 0; j + 1 < n; ++j) {
        if (std::abs(diag[j]) >= std::abs(lower[j])) {
            if (diag[j] == C{}) throw std::domain_error("singular Helmholtz matrix");
            const C f = lower[j] / diag[j];
            diag[j+1] -= f * upper[j];
            rhs[j+1] -= f * rhs[j];
            lower[j] = {};
        } else {
            const C f = diag[j] / lower[j];
            diag[j] = lower[j];
            const C next = diag[j+1];
            diag[j+1] = upper[j] - f * next;
            upper[j] = next;
            lower[j] = {};
            if (j + 2 < n) {
                lower[j] = upper[j+1];
                upper[j+1] = -f * lower[j];
            }
            std::swap(rhs[j], rhs[j+1]);
            rhs[j+1] -= f * rhs[j];
        }
        check(diag[j+1]); check(rhs[j+1]);
    }
    for (std::size_t j = n; j-- > 0;) {
        if (diag[j] == C{}) throw std::domain_error("singular Helmholtz matrix");
        if (j + 1 < n) rhs[j] -= upper[j] * rhs[j+1];
        if (j + 2 < n) rhs[j] -= lower[j] * rhs[j+2];
        rhs[j] /= diag[j];
        check(rhs[j]);
    }
    return rhs;
}
} // namespace

Helmholtz1DResult solve_helmholtz_1d(double wavelength, double ni, double nt,
    std::span<const double> nodes, std::span<const std::complex<double>> eps) {
    positive(wavelength); positive(ni); positive(nt);
    if (nodes.size() < 2 || eps.size() != nodes.size()-1)
        throw std::invalid_argument("inconsistent Helmholtz mesh and materials");
    const double k = 2.0 * std::numbers::pi_v<double> / wavelength;
    positive(k);
    const auto n = nodes.size();
    std::vector<C> diag(n), off(n-1), rhs(n);
    std::vector<double> h(n-1);
    for (std::size_t j=0; j+1<n; ++j) {
        if (!std::isfinite(nodes[j]) || !std::isfinite(nodes[j+1]) ||
            !finite(eps[j]) || eps[j].imag() < 0.0)
            throw std::invalid_argument("invalid mesh coordinate or passive permittivity");
        // Dimensionless coordinate s=k0*x: avoids SI stiffness scaling.
        h[j] = k * (nodes[j+1] - nodes[j]);
        positive(h[j]);
        const C mass = eps[j] * (h[j]/6.0);
        const C diagonal = 1.0/h[j] - 2.0*mass;
        off[j] = -1.0/h[j] - mass;
        diag[j] += diagonal; diag[j+1] += diagonal;
        check(off[j]); check(diag[j]); check(diag[j+1]);
    }
    diag.front() -= I*ni; diag.back() -= I*nt;
    rhs.front() = -2.0*I*ni;
    check(diag.front()); check(diag.back()); check(rhs.front());
    auto field = solve(off, diag, off, std::move(rhs));
    const C r = field.front()-1.0, t = field.back();
    const double R = std::norm(r), T = (nt/ni)*std::norm(t);
    std::vector<double> loss(n-1);
    double A=0.0;
    for (std::size_t j=0; j+1<n; ++j) {
        // Exact integral of |E_h|^2 for the piecewise-linear field.
        const C a=field[j], b=field[j+1];
        const double integral = h[j]/3.0 *
            (std::norm(a)+std::norm(b)+std::real(std::conj(a)*b));
        loss[j] = eps[j].imag()/ni * integral;
        A += loss[j];
    }
    const double residual = 1.0-R-T-A;
    if (!std::isfinite(R) || !std::isfinite(T) ||
        !std::isfinite(A) || !std::isfinite(residual))
        throw std::overflow_error("nonfinite Helmholtz power balance");
    return {std::move(field), std::move(loss), r, t, R, T, A, residual};
}

ObliqueHelmholtz1DResult solve_oblique_helmholtz_1d(
    double wavelength, double ni, double nt, double angle, Polarization pol,
    std::span<const double> nodes, std::span<const C> eps) {
    positive(wavelength); positive(ni); positive(nt);
    if (!std::isfinite(angle) || angle < 0.0 ||
        angle >= std::numbers::pi_v<double>/2 ||
        (pol != Polarization::te && pol != Polarization::tm) ||
        nodes.size()<2 || eps.size()!=nodes.size()-1)
        throw std::invalid_argument("invalid oblique Helmholtz inputs");
    const bool tm=pol==Polarization::tm;
    const double k=2.0*std::numbers::pi_v<double>/wavelength;
    positive(k);
    const double beta=ni*std::sin(angle), b2=beta*beta;
    const double ei=ni*ni, et=nt*nt;
    positive(ei); positive(et);
    if (!std::isfinite(b2)) throw std::overflow_error("nonfinite transverse wavevector");
    const double qi=ni*std::cos(angle);
    const C qt=std::sqrt(C{et-b2,0.0});
    if (qt==C{}) throw std::domain_error("exact exit critical angle unsupported");
    // Scalar flux coefficient: q for TE and q/epsilon for TM.
    const double ai=tm?qi/ei:qi;
    const C at=tm?qt/et:qt;
    positive(ai); check(at);
    const auto n=nodes.size();
    std::vector<C> diag(n), off(n-1), rhs(n);
    std::vector<double> h(n-1), loss_factor(n-1);
    for (std::size_t j=0;j+1<n;++j) {
        if (!std::isfinite(nodes[j]) || !std::isfinite(nodes[j+1]) ||
            !finite(eps[j]) || eps[j].imag()<0.0 || (tm && eps[j]==C{}))
            throw std::invalid_argument("invalid oblique mesh or passive permittivity");
        h[j]=k*(nodes[j+1]-nodes[j]);
        positive(h[j]);
        const C a=tm?1.0/eps[j]:C{1.0,0.0};
        const C c=a*(eps[j]-b2);
        check(a); check(c);
        loss_factor[j]=tm?-a.imag():eps[j].imag();
        const C mass=c*(h[j]/6.0);
        const C d=a/h[j]-2.0*mass;
        off[j]=-a/h[j]-mass;
        diag[j]+=d; diag[j+1]+=d;
        check(off[j]); check(diag[j]); check(diag[j+1]);
    }
    diag.front()-=I*ai;
    diag.back()-=I*at;
    rhs.front()=-2.0*I*ai;
    check(diag.front()); check(diag.back()); check(rhs.front());
    auto u=solve(off,diag,off,std::move(rhs));
    const C rs=u.front()-1.0, ts=u.back();
    // The backward magnetic wave reverses the tangential electric sign.
    const C r=tm?-rs:rs, t=tm?ts*(at/ai):ts;
    check(r); check(t);
    const double R=std::norm(rs), T=(at.real()/ai)*std::norm(ts);
    std::vector<double> loss(n-1);
    double A=0;
    for (std::size_t j=0;j+1<n;++j) {
        const C a=u[j], b=u[j+1];
        const double integral=h[j]/3.0*
            (std::norm(a)+std::norm(b)+std::real(std::conj(a)*b));
        // TM dissipates through both transverse and longitudinal E.
        const double electric_integral=tm?
            std::norm(b-a)/h[j]+b2*integral:integral;
        loss[j]=loss_factor[j]/ai*electric_integral;
        if (!std::isfinite(loss[j]) || loss[j]<0.0)
            throw std::overflow_error("invalid oblique cell absorption");
        A+=loss[j];
    }
    const double residual=1.0-R-T-A;
    if (!std::isfinite(R) || !std::isfinite(T) ||
        !std::isfinite(A) || !std::isfinite(residual))
        throw std::overflow_error("nonfinite oblique power balance");
    return {std::move(u),std::move(loss),r,t,R,T,A,residual};
}
} // namespace nova::optics
