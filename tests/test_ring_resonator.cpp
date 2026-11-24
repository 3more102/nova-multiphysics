#include "nova/physics/ring_resonator.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace {
void near(double x, double y, double tolerance, const char* message) {
    if (!std::isfinite(x) || std::abs(x - y) > tolerance)
        throw std::runtime_error(message);
}
}

int main() {
    using nova::optics::all_pass_ring;
    constexpr double pi = std::numbers::pi_v<double>;
    try {
        const auto critical = all_pass_ring(0.9, 0.9, 0.0);
        near(critical.through_power, 0.0, 1e-14, "critical coupling");
        near(critical.dissipated_power, 1.0, 1e-14, "critical dissipation");
        near(critical.circulating_power, 1.0 / (1.0 - 0.81), 1e-12, "critical buildup");
        const auto anti = all_pass_ring(0.9, 0.8, pi);
        near(anti.through_power, std::pow(1.7 / 1.72, 2), 1e-14, "antiresonance");
        near(all_pass_ring(0.0, 0.7, 0.2).through_power, 0.49, 1e-14, "full coupling");
        near(all_pass_ring(1.0, 1.0, 0.0).circulating_power, 0.0, 0.0, "disconnected");
        near(all_pass_ring(0.8, 0.0, 0.3).through_power, 0.64, 1e-14, "zero survival");
        near(all_pass_ring(0.9, 1.0, 0.0).through_power, 1.0, 1e-14, "lossless");
        const double almost_one = std::nextafter(1.0, 0.0);
        const auto narrow = all_pass_ring(almost_one, almost_one, 0.0);
        near(narrow.through_power, 0.0, 1e-14, "near-unit critical coupling");
        near(narrow.dissipated_power, 1.0, 1e-14, "near-unit loss");
        for (int i = 0; i <= 10; ++i) {
            for (int j = 0; j <= 10; ++j) {
                const double t = i / 10.0, a = j / 10.0;
                for (int k = -20; k <= 20; ++k) {
                    const double p = k * pi / 20.0;
                    const auto r = all_pass_ring(t, a, p);
                    near(r.power_balance_defect, 0.0, 2e-14, "passive power balance");
                    near(r.through_power, all_pass_ring(t, a, p + 2.0*pi).through_power,
                         2e-14, "periodicity");
                    if (r.through_power < 0.0 || r.dissipated_power < 0.0)
                        throw std::runtime_error("negative power");
                }
            }
        }
        int rejected = 0;
        for (const double invalid : {-0.1, 1.1, std::numeric_limits<double>::infinity(),
                                     std::numeric_limits<double>::quiet_NaN()}) {
            try { (void)all_pass_ring(invalid, 0.8, 0.0); }
            catch (const std::invalid_argument&) { ++rejected; }
            try { (void)all_pass_ring(0.8, invalid, 0.0); }
            catch (const std::invalid_argument&) { ++rejected; }
        }
        try { (void)all_pass_ring(0.8, 0.8, std::numeric_limits<double>::infinity()); }
        catch (const std::invalid_argument&) { ++rejected; }
        if (rejected != 9) throw std::runtime_error("invalid inputs accepted");
        std::cout << "NOVA ring reference: 4961 parameter combinations passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
