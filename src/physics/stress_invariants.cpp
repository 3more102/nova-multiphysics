#include "nova/physics/linear_elasticity_3d.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace nova::physics {
namespace {
double normalize(ElasticTensor6& s) {
    double magnitude=0;
    for (double x:s) {
        if (!std::isfinite(x)) throw std::invalid_argument("Stress tensor must be finite.");
        magnitude=std::max(magnitude,std::abs(x));
    }
    if (magnitude) for (double& x:s) x/=magnitude;
    return magnitude;
}
}
double von_mises_stress_3d(const ElasticTensor6& stress) {
    auto s=stress;
    const double magnitude=normalize(s);
    const double a=s[0]-s[1], b=s[1]-s[2], c=s[2]-s[0];
    return magnitude*std::sqrt(.5*(a*a+b*b+c*c)+3*(s[3]*s[3]+s[4]*s[4]+s[5]*s[5]));
}
ElasticVector3D principal_stresses_3d(const ElasticTensor6& stress) {
    auto s=stress;
    const double magnitude=normalize(s);
    if (s[3]==0 && s[4]==0 && s[5]==0) {
        ElasticVector3D result{stress[0],stress[1],stress[2]};
        std::sort(result.rbegin(),result.rend());
        return result;
    }
    const double mean=(s[0]+s[1]+s[2])/3;
    for (unsigned i=0;i<3;++i) s[i]-=mean;
    const double p=std::hypot(std::hypot(s[0],s[1],s[2]),std::sqrt(2.)*std::hypot(s[3],s[4],s[5]))/std::sqrt(6.);
    if (p==0) return {mean*magnitude,mean*magnitude,mean*magnitude};
    for (double& x:s) x/=p;
    const double det=s[0]*(s[1]*s[2]-s[4]*s[4])-s[3]*(s[3]*s[2]-s[4]*s[5])+s[5]*(s[3]*s[4]-s[1]*s[5]);
    const double angle=std::acos(std::clamp(det/2,-1.,1.))/3;
    const double high=2*p*std::cos(angle), low=2*p*std::cos(angle+2*std::numbers::pi/3);
    return {magnitude*(mean+high),magnitude*(mean-high-low),magnitude*(mean+low)};
}
}
