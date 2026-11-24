#include "nova/physics/ring_resonator.hpp"
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {
void near(double actual,double expected,double tol=1e-12) {
    if (!std::isfinite(actual) || std::abs(actual-expected)>tol)
        throw std::runtime_error("ring reference mismatch");
}
}
int main() {
    using nova::optics::add_drop_ring;
    try {
        auto r=add_drop_ring(0.8,0.8,1.0,0.0);
        near(r.through_power,0); near(r.drop_power,1); near(r.absorbed_power,0);
        near(r.circulating_power,1.0/(1.0-0.64));
        r=add_drop_ring(0.8,1.0,0.8,0.0);
        near(r.through_power,0); near(r.drop_power,0); near(r.absorbed_power,1);
        r=add_drop_ring(1.0,1.0,1.0,0.0);
        near(r.through_power,1); near(r.circulating_power,0);
        r=add_drop_ring(0.0,0.0,0.7,0.4);
        near(r.through_power,0); near(r.drop_power,0.7); near(r.absorbed_power,0.3);
        r=add_drop_ring(0.99999999,0.99999999,1.0,0.0);
        near(r.drop_power,1,1e-8);
        for(double a:{0.0,0.2,0.9,1.0})
            for(double t1:{0.0,0.3,0.95,1.0})
                for(double t2:{0.0,0.6,0.99,1.0})
                    for(int i=0;i<33;++i) {
                        const double p=i*std::numbers::pi_v<double>/16;
                        r=add_drop_ring(t1,t2,a,p);
                        near(r.through_power+r.drop_power+r.absorbed_power,1);
                        const auto shifted=add_drop_ring(t1,t2,a,p+2*std::numbers::pi_v<double>);
                        near(r.through_power,shifted.through_power);
                        near(r.drop_power,shifted.drop_power);
                    }
        bool rejected=false;
        try { (void)add_drop_ring(1.1,0.5,1,0); }
        catch(const std::invalid_argument&) { rejected=true; }
        if(!rejected) throw std::runtime_error("invalid coupling accepted");
        std::cout<<"Ring reference checks passed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
