#include "nova/physics/microring.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace {
void near(double x, double y, double tol, const char* message) {
    if (!std::isfinite(x) || std::abs(x-y)>tol) throw std::runtime_error(message);
}
}
int main() {
    using nova::optics::microring_response;
    try {
        const double pi = std::numbers::pi_v<double>;
        const auto critical = microring_response(0.9,1.0,0.9,0.0);
        near(critical.through_power,0.0,1e-14,"critical coupling extinction");
        near(critical.absorbed_power,1.0,1e-13,"critical absorption");
        const auto matched = microring_response(0.8,0.8,1.0,0.0);
        near(matched.through_power,0.0,1e-14,"matched through extinction");
        near(matched.drop_power,1.0,1e-13,"matched drop transfer");
        const auto anti = microring_response(0.8,0.8,1.0,pi);
        near(anti.through_power,std::pow(1.6/1.64,2),1e-13,"antiresonance");
        const auto uncoupled = microring_response(1.0,1.0,1.0,0.0);
        near(uncoupled.through_power,1.0,0.0,"uncoupled singular limit");
        const auto opaque = microring_response(0.6,0.7,0.0,0.3);
        near(opaque.through_power,0.36,1e-14,"opaque ring through");
        near(opaque.absorbed_power,0.64,1e-14,"opaque ring absorption");
        const double almost = std::nextafter(1.0,0.0);
        const auto high = microring_response(almost,1.0,almost,0.0);
        near(high.through_power,0.0,1e-14,"high finesse critical coupling");
        near(high.absorbed_power,1.0,1e-12,"high finesse absorption");
        for (double t1 : {0.0,0.2,0.8,0.99,1.0})
            for (double t2 : {0.0,0.3,0.9,1.0})
                for (double a : {0.0,0.5,0.999,1.0})
                    for (int i=0;i<=32;++i) {
                        const double phase = i*pi/16.0;
                        const auto x = microring_response(t1,t2,a,phase);
                        near(x.through_power+x.drop_power+x.absorbed_power,
                             1.0,1e-12,"passive power conservation");
                        const auto periodic=microring_response(t1,t2,a,phase+2.0*pi);
                        near(x.through_power,periodic.through_power,1e-12,"periodic through");
                        near(x.drop_power,periodic.drop_power,1e-12,"periodic drop");
                        if (a==1.0) near(x.absorbed_power,0.0,0.0,"lossless absorption");
                    }
        for (double invalid : {-0.1,1.1,std::numeric_limits<double>::quiet_NaN()}) {
            bool rejected=false;
            try { (void)microring_response(invalid,0.5,0.9,0.0); }
            catch (const std::invalid_argument&) { rejected=true; }
            if (!rejected) throw std::runtime_error("invalid coupling accepted");
        }
        std::cout << "NOVA microring reference cases passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
