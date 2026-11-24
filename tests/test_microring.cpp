#include "nova/physics/microring.hpp"
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {
void near(double x,double y,double tol) {
    if (!std::isfinite(x) || std::abs(x-y)>tol)
        throw std::runtime_error("microring reference mismatch");
}
}
int main() {
    using nova::optics::microring;
    try {
        auto critical = microring(0.9,1,0.9,0);
        near(critical.through_power,0,1e-14);
        near(critical.absorbed_power,1,1e-13);
        auto add_drop = microring(0.8,0.8,1,0);
        near(add_drop.through_power,0,1e-14);
        near(add_drop.drop_power,1,1e-13);
        near(std::abs(microring(0.8,0.8,1,2*std::numbers::pi_v<double>).drop_amplitude
            + add_drop.drop_amplitude),0,1e-12);
        auto uncoupled = microring(1,0.8,0.9,0.3);
        near(uncoupled.through_power,1,1e-14);
        near(uncoupled.circulating_power,0,1e-14);
        auto extinguished = microring(0.6,0.7,0,0.3);
        near(extinguished.through_power,0.36,1e-14);
        near(extinguished.absorbed_power,0.64,1e-14);
        const double t=0.8, a=0.9;
        for (int i=0;i<=1000;++i) {
            const double p=2*std::numbers::pi_v<double>*i/1000;
            const auto r=microring(t,1,a,p);
            near(r.through_power,(t*t+a*a-2*t*a*std::cos(p))/
                (1+t*t*a*a-2*t*a*std::cos(p)),1e-12);
            near(r.power_balance_defect,0,1e-12);
            near(microring(t,1,1,p).through_power,1,1e-12);
            near(microring(0.8,0.7,0.93,p).power_balance_defect,0,1e-12);
        }
        for (double p : {0.0, 2*std::numbers::pi_v<double>}) {
            bool rejected=false;
            try { (void)microring(1,1,1,p); }
            catch(const std::domain_error&) { rejected=true; }
            if (!rejected) return 1;
        }
        bool rejected=false;
        try { (void)microring(1.1,1,1,0); }
        catch(const std::invalid_argument&) { rejected=true; }
        if (!rejected) return 1;
        std::cout << "Microring reference tests passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
