#include "nova/physics/optics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
using nova::optics::Complex;
using nova::optics::Layer;
using nova::optics::StackResult;
using Matrix = std::array<Complex, 4>;
constexpr Complex i{0.0, 1.0};
constexpr double tol = 2e-11;

Matrix multiply(const Matrix& a, const Matrix& b) {
    return {a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3]};
}
// Independent E,H characteristic matrix, mapping right boundary to left.
// Used only for moderate optical depth; unlike production scattering recursion,
// this reference is not suitable for extremely opaque stacks.
std::array<Complex,2> reference(double wavelength, double n0, double ns,
                                 const std::vector<Layer>& layers) {
    Matrix m{1.0, 0.0, 0.0, 1.0};
    for (const auto& l : layers) {
        const auto delta = 2.0*std::numbers::pi_v<double>*l.index*l.thickness_m/wavelength;
        const auto c = std::cos(delta), s = std::sin(delta);
        m = multiply(m, {c, -i*s/l.index, -i*l.index*s, c});
    }
    const auto b = m[0]+ns*m[1], h = m[2]+ns*m[3];
    return {(n0*b-h)/(n0*b+h), 2.0*n0/(n0*b+h)};
}
void near(Complex a, Complex b, const char* message) {
    if (!(std::abs(a-b) <= tol)) throw std::runtime_error(message);
}
void check(const StackResult& a, const StackResult& b) {
    near(a.reflection,b.reflection,"reflection invariance");
    near(a.transmission,b.transmission,"transmission invariance");
}
}
int main() {
    using nova::optics::normal_incidence_stack;
    try {
        std::size_t cases = 0;
        for (double wavelength : {0.45e-6, 0.65e-6, 1.31e-6, 1.55e-6}) {
            for (double n0 : {1.0, 1.33, 2.0}) {
                for (double ns : {1.0, 1.52, 2.4}) {
                    for (double loss : {0.0, 0.03, 0.2}) {
                        const std::vector<Layer> layers{{{1.7,loss},110e-9},
                            {{2.2,0.5*loss},90e-9},{{1.2,loss},230e-9}};
                        const auto a = normal_incidence_stack(wavelength,n0,ns,layers);
                        const auto ref = reference(wavelength,n0,ns,layers);
                        near(a.reflection,ref[0],"characteristic matrix reflection");
                        near(a.transmission,ref[1],"characteristic matrix transmission");
                        if (!(a.absorptance >= -tol)) throw std::runtime_error("passivity");
                        if (loss == 0.0) near(a.reflectance+a.transmittance,1.0,"lossless conservation");
                        auto reversed = layers;
                        std::reverse(reversed.begin(),reversed.end());
                        const auto rev = normal_incidence_stack(wavelength,ns,n0,reversed);
                        near(a.transmittance,rev.transmittance,"reciprocal power transmission");
                        if (loss == 0.0) near(a.reflectance,rev.reflectance,"lossless reciprocal reflection");
                        auto split = layers;
                        split[0].thickness_m *= 0.5;
                        split.insert(split.begin(),split.front());
                        check(a,normal_incidence_stack(wavelength,n0,ns,split));
                        auto zero = layers;
                        zero.insert(zero.begin()+1,Layer{{3.1,0.4},0.0});
                        check(a,normal_incidence_stack(wavelength,n0,ns,zero));
                        ++cases;
                    }
                }
            }
        }
        std::cout << "NOVA independent optics crosscheck: " << cases << " parameter cases passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
