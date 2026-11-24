#include "nova/linalg/incomplete_lu.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
namespace nova::linalg {
IncompleteLU0::IncompleteLU0(const CsrMatrix& matrix) {
    const auto n=matrix.rows();
    if (n==0 || n!=matrix.columns()) throw std::invalid_argument("ILU(0) requires a nonempty square system.");
    offsets_.resize(n+1);
    diagonal_.assign(n,matrix.non_zeros());
    columns_.reserve(matrix.non_zeros()); values_.reserve(matrix.non_zeros());
    matrix.for_each_nonzero([&](auto i,auto j,double x){
        if (i==j) diagonal_[i]=values_.size();
        ++offsets_[i+1]; columns_.push_back(j); values_.push_back(x);
    });
    for (std::size_t i=0;i<n;++i) offsets_[i+1]+=offsets_[i];
    std::vector<std::size_t> positions(n,matrix.non_zeros());
    for (std::size_t i=0;i<n;++i) {
        if (diagonal_[i]==matrix.non_zeros()) throw std::invalid_argument("ILU(0) missing diagonal.");
        double magnitude=0;
        for (auto k=offsets_[i];k<offsets_[i+1];++k) {
            positions[columns_[k]]=k;
            magnitude=std::max(magnitude,std::abs(values_[k]));
        }
        for (auto k=offsets_[i];k<diagonal_[i];++k) {
            const auto j=columns_[k];
            const double l=values_[k]/values_[diagonal_[j]];
            values_[k]=l;
            for (auto p=diagonal_[j]+1;p<offsets_[j+1];++p) {
                const auto dest=positions[columns_[p]];
                if (dest!=matrix.non_zeros()) values_[dest]-=l*values_[p];
            }
        }
        if (!(std::abs(values_[diagonal_[i]])>32*std::numeric_limits<double>::epsilon()*magnitude))
            throw std::domain_error("ILU(0) unresolved pivot at row "+std::to_string(i));
        for (auto k=offsets_[i];k<offsets_[i+1];++k) {
            if (!std::isfinite(values_[k])) throw std::domain_error("Non-finite ILU(0) factor.");
            positions[columns_[k]]=matrix.non_zeros();
        }
    }
}
std::vector<double> IncompleteLU0::solve(std::span<const double> rhs) const {
    if (rhs.size()!=diagonal_.size()) throw std::invalid_argument("ILU right-hand side dimension mismatch.");
    std::vector<double> x(rhs.begin(),rhs.end());
    for (std::size_t i=0;i<x.size();++i)
        for (auto k=offsets_[i];k<diagonal_[i];++k) x[i]-=values_[k]*x[columns_[k]];
    for (std::size_t i=x.size();i-->0;) {
        for (auto k=diagonal_[i]+1;k<offsets_[i+1];++k) x[i]-=values_[k]*x[columns_[k]];
        x[i]/=values_[diagonal_[i]];
    }
    return x;
}
}
