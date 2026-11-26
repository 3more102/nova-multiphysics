#include "nova/linalg/symmetric_eigen.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace nova::linalg {
namespace {
struct Dense {
    std::size_t n;
    std::vector<double> data;
    explicit Dense(std::size_t size):n(size),data(size*size){}
    double& operator()(std::size_t i,std::size_t j){return data[i*n+j];}
    double operator()(std::size_t i,std::size_t j) const{return data[i*n+j];}
};
double norm(const std::vector<double>& x) {
    double out=0.0;
    for (double value:x) out=std::hypot(out,value);
    return out;
}
// Nonnegative magnitudes represented as mantissa * 2^exponent. Keep
// residual normalization outside the double exponent range until the
// final ratio, including when the Frobenius norm itself is too large.
struct Magnitude {
    double mantissa{};
    int exponent{};
};
Magnitude magnitude(double value) {
    int exponent=0;
    const double mantissa=std::frexp(value,&exponent);
    return {mantissa,exponent};
}
Magnitude multiply(Magnitude a,Magnitude b) {
    if (a.mantissa==0.0 || b.mantissa==0.0) return {};
    auto result=magnitude(a.mantissa*b.mantissa);
    result.exponent+=a.exponent+b.exponent;
    return result;
}
Magnitude add(Magnitude a,Magnitude b) {
    if (a.mantissa==0.0) return b;
    if (b.mantissa==0.0) return a;
    if (a.exponent<b.exponent) std::swap(a,b);
    auto result=magnitude(a.mantissa+std::scalbn(b.mantissa,b.exponent-a.exponent));
    result.exponent+=a.exponent;
    return result;
}
Magnitude scaled_norm(const std::vector<double>& values) {
    double scale=0.0;
    for (double value:values) {
        if (!std::isfinite(value)) throw std::overflow_error("Eigen residual norm overflow.");
        scale=std::max(scale,std::abs(value));
    }
    if (scale==0.0) return {};
    double sum=0.0;
    for (double value:values) {
        const double normalized=value/scale;
        sum+=normalized*normalized;
    }
    return multiply(magnitude(scale),magnitude(std::sqrt(sum)));
}
double relative_residual(const std::vector<double>& r,const std::vector<double>& x,
                         double lambda,Magnitude kn,Magnitude mn) {
    const auto numerator=scaled_norm(r);
    const auto denominator=multiply(add(kn,multiply(magnitude(std::abs(lambda)),mn)),
                                    scaled_norm(x));
    if (denominator.mantissa==0.0) return std::scalbn(numerator.mantissa,numerator.exponent);
    if (numerator.mantissa==0.0) return 0.0;
    return std::scalbn(numerator.mantissa/denominator.mantissa,
                       numerator.exponent-denominator.exponent);
}
Dense dense_symmetric(const CsrMatrix& a) {
    Dense out(a.rows());
    a.for_each_nonzero([&](std::size_t i,std::size_t j,double v){out(i,j)=v;});
    double scale=0.0;
    for (double v:out.data) {
        if (!std::isfinite(v)) throw std::invalid_argument("Nonfinite eigen matrix.");
        scale=std::max(scale,std::abs(v));
    }
    for (std::size_t i=0;i<out.n;++i)
        for (std::size_t j=i+1;j<out.n;++j) {
            if (std::abs(out(i,j)-out(j,i))>1e-13*scale)
                throw std::invalid_argument("Eigen matrix must be symmetric.");
            const double v=0.5*out(i,j)+0.5*out(j,i);
            out(i,j)=out(j,i)=v;
        }
    return out;
}
Dense cholesky(const Dense& m) {
    Dense l(m.n);
    for (std::size_t i=0;i<m.n;++i)
        for (std::size_t j=0;j<=i;++j) {
            double v=m(i,j);
            for (std::size_t k=0;k<j;++k) v-=l(i,k)*l(j,k);
            if (i==j) {
                if (!(v>0.0) || !std::isfinite(v))
                    throw std::invalid_argument("Eigen mass is not numerically SPD.");
                l(i,j)=std::sqrt(v);
            } else l(i,j)=v/l(j,j);
            if (!std::isfinite(l(i,j))) throw std::overflow_error("Cholesky overflow.");
        }
    return l;
}
Dense whiten(const Dense& k,const Dense& l) {
    Dense x(k.n),b(k.n);
    for (std::size_t j=0;j<k.n;++j)
        for (std::size_t i=0;i<k.n;++i) {
            double v=k(i,j);
            for (std::size_t p=0;p<i;++p) v-=l(i,p)*x(p,j);
            x(i,j)=v/l(i,i);
        }
    for (std::size_t i=0;i<k.n;++i)
        for (std::size_t j=0;j<k.n;++j) {
            double v=x(i,j);
            for (std::size_t p=0;p<j;++p) v-=l(j,p)*b(i,p);
            b(i,j)=v/l(j,j);
        }
    for (std::size_t i=0;i<k.n;++i)
        for (std::size_t j=i;j<k.n;++j) {
            const double v=0.5*b(i,j)+0.5*b(j,i);
            if (!std::isfinite(v)) throw std::overflow_error("Eigen reduction overflow.");
            b(i,j)=b(j,i)=v;
        }
    return b;
}
double off_norm(const Dense& b) {
    double result=0.0;
    for (std::size_t i=0;i<b.n;++i)
        for (std::size_t j=i+1;j<b.n;++j) result=std::hypot(result,b(i,j));
    return result;
}
}

SymmetricEigenResult solve_symmetric_generalized_eigen(
    const CsrMatrix& stiffness,const CsrMatrix& mass,const SymmetricEigenOptions& options) {
    const auto n=stiffness.rows();
    if (n==0 || n!=stiffness.columns() || n!=mass.rows() || n!=mass.columns()
        || n>options.max_dimension || n>512 || options.max_sweeps==0
        || !(options.relative_tolerance>0.0) || options.relative_tolerance>=1.0
        || !std::isfinite(options.relative_tolerance))
        throw std::invalid_argument("Invalid generalized eigen dimensions or options.");
    const auto k=dense_symmetric(stiffness),m=dense_symmetric(mass);
    const auto l=cholesky(m);
    auto b=whiten(k,l);
    Dense v(n);
    for (std::size_t i=0;i<n;++i) v(i,i)=1.0;
    double scale=0.0;
    for (double value:b.data) scale=std::max(scale,std::abs(value));
    if (scale>0.0) for (double& value:b.data) value/=scale;
    const double goal=options.relative_tolerance*norm(b.data);
    SymmetricEigenResult result;
    for (;result.sweeps<options.max_sweeps && off_norm(b)>goal;++result.sweeps) {
        const double skip=goal/static_cast<double>(n);
        for (std::size_t p=0;p<n;++p)
            for (std::size_t q=p+1;q<n;++q) {
                const double apq=b(p,q);
                if (std::abs(apq)<=skip) continue;
                const double tau=(b(q,q)-b(p,p))/(2.0*apq);
                const double t=std::copysign(1.0,tau)/(std::abs(tau)+std::hypot(1.0,tau));
                const double c=1.0/std::hypot(1.0,t),s=t*c;
                b(p,p)-=t*apq; b(q,q)+=t*apq;
                b(p,q)=b(q,p)=0.0;
                for (std::size_t i=0;i<n;++i) {
                    if (i!=p && i!=q) {
                        const double a=b(i,p),z=b(i,q);
                        b(i,p)=b(p,i)=c*a-s*z;
                        b(i,q)=b(q,i)=s*a+c*z;
                    }
                    const double a=v(i,p),z=v(i,q);
                    v(i,p)=c*a-s*z; v(i,q)=s*a+c*z;
                }
            }
    }
    result.converged=off_norm(b)<=goal;
    const auto kn=scaled_norm(k.data),mn=scaled_norm(m.data);
    for (std::size_t j=0;j<n;++j) {
        SymmetricEigenpair pair;
        pair.eigenvalue=b(j,j)*scale;
        auto& x=pair.eigenvector;
        x.resize(n);
        for (std::size_t i=n;i-->0;) {
            double z=v(i,j);
            for (std::size_t p=i+1;p<n;++p) z-=l(p,i)*x[p];
            x[i]=z/l(i,i);
        }
        const auto mx=mass.multiply(x);
        double mass_norm=0.0;
        for (std::size_t i=0;i<n;++i) mass_norm+=x[i]*mx[i];
        if (!(mass_norm>0.0)||!std::isfinite(mass_norm)||!std::isfinite(pair.eigenvalue))
            throw std::overflow_error("Invalid eigenpair normalization.");
        const double normalization=std::sqrt(mass_norm);
        for (double& value:x) value/=normalization;
        auto r=stiffness.multiply(x);
        const auto mx_normalized=mass.multiply(x);
        for (std::size_t i=0;i<n;++i) r[i]-=pair.eigenvalue*mx_normalized[i];
        pair.relative_residual=relative_residual(r,x,pair.eigenvalue,kn,mn);
        if (!std::isfinite(pair.relative_residual))
            throw std::overflow_error("Eigen residual overflow.");
        if (pair.relative_residual>std::max(1e-10,10*options.relative_tolerance))
            result.converged=false;
        result.eigenpairs.push_back(std::move(pair));
    }
    std::sort(result.eigenpairs.begin(),result.eigenpairs.end(),
        [](const auto& a,const auto& z){return a.eigenvalue<z.eigenvalue;});
    return result;
}

} // namespace nova::linalg
