#include "test_framework.hpp"
#include "nova/physics/floating_electrodes_3d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using nova::physics::CapacitanceMatrix3D;
using nova::physics::FixedElectrodePotential3D;
using nova::physics::FloatingElectrodeCharge3D;
using nova::physics::solve_floating_electrode_potentials_3d;

template <class Callable>
void require_invalid(const Callable& operation) {
    bool invalid = false;
    try { operation(); }
    catch (const std::invalid_argument&) { invalid = true; }
    NOVA_REQUIRE(invalid);
}

} // namespace

NOVA_TEST(floating_electrodes_parallel_plate_charge_to_voltage) {
    const CapacitanceMatrix3D c{{10, 20}, {2.0, -2.0, -2.0, 2.0}};
    const std::vector<FixedElectrodePotential3D> fixed{{10, 0.0}};
    const std::vector<FloatingElectrodeCharge3D> free{{20, 20.0}};
    const auto result = solve_floating_electrode_potentials_3d(c, fixed, free);
    NOVA_REQUIRE(result.electrode_ids == c.electrode_ids);
    NOVA_REQUIRE(result.potential_volts.size() == 2U);
    NOVA_REQUIRE_NEAR(result.potential_volts[0], 0.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.potential_volts[1], 10.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[0], -20.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[1], 20.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.max_floating_charge_residual_coulombs, 0.0, 1.0e-12);

    const auto shifted = solve_floating_electrode_potentials_3d(
        c, std::vector<FixedElectrodePotential3D>{{20, 4.0}},
        std::vector<FloatingElectrodeCharge3D>{{10, -6.0}});
    NOVA_REQUIRE_NEAR(shifted.potential_volts[0], 1.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(shifted.potential_volts[1], 4.0, 1.0e-12);
}

NOVA_TEST(floating_electrodes_three_terminal_coupling) {
    const CapacitanceMatrix3D c{
        {10, 20, 30},
        {3.0, -1.0, -2.0, -1.0, 2.0, -1.0, -2.0, -1.0, 3.0}
    };
    // Deliberately supply floating terminals in reverse electrode-ID order.
    const auto result = solve_floating_electrode_potentials_3d(
        c, std::vector<FixedElectrodePotential3D>{{10, 3.0}},
        std::vector<FloatingElectrodeCharge3D>{{30, -4.0}, {20, 3.0}});
    NOVA_REQUIRE_NEAR(result.potential_volts[0], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.potential_volts[1], 4.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.potential_volts[2], 2.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[0], 1.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[1], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[2], -4.0, 1.0e-12);
}

NOVA_TEST(floating_electrodes_external_ground_and_all_fixed) {
    // A positive-definite terminal matrix models an independent grounded
    // reference; unlike an isolated Maxwell matrix, all terminals may float.
    const CapacitanceMatrix3D c{{1, 2}, {3.0, -1.0, -1.0, 2.0}};
    const auto all_free = solve_floating_electrode_potentials_3d(
        c, std::vector<FixedElectrodePotential3D>{},
        std::vector<FloatingElectrodeCharge3D>{{1, 3.0}, {2, 4.0}});
    NOVA_REQUIRE_NEAR(all_free.potential_volts[0], 2.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_free.potential_volts[1], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_free.estimated_charge_coulombs[0], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_free.estimated_charge_coulombs[1], 4.0, 1.0e-12);

    const auto all_fixed = solve_floating_electrode_potentials_3d(
        c, std::vector<FixedElectrodePotential3D>{{2, 3.0}, {1, 2.0}},
        std::vector<FloatingElectrodeCharge3D>{});
    NOVA_REQUIRE_NEAR(all_fixed.estimated_charge_coulombs[0], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_fixed.estimated_charge_coulombs[1], 4.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_fixed.max_floating_charge_residual_coulombs, 0.0, 1.0e-12);
}

NOVA_TEST(floating_electrodes_rejects_invalid_inputs) {
    const CapacitanceMatrix3D c{{10, 20}, {2.0, -2.0, -2.0, 2.0}};
    const std::vector<FixedElectrodePotential3D> fixed{{10, 0.0}};
    const std::vector<FloatingElectrodeCharge3D> free{{20, 20.0}};
    const auto run = [&](const CapacitanceMatrix3D& matrix,
                         const std::vector<FixedElectrodePotential3D>& k,
                         const std::vector<FloatingElectrodeCharge3D>& f) {
        (void)solve_floating_electrode_potentials_3d(matrix, k, f);
    };

    require_invalid([&] { run({{}, {}}, fixed, free); });
    require_invalid([&] { run({{10, 20}, {2.0}}, fixed, free); });
    require_invalid([&] { run({{20, 10}, c.matrix_farads}, fixed, free); });
    require_invalid([&] { run({{10, 10}, c.matrix_farads}, fixed, free); });
    require_invalid([&] { run({{-1, 20}, c.matrix_farads}, fixed, free); });
    require_invalid([&] {
        run({{10, 20}, {2.0, -2.0, -2.0, std::numeric_limits<double>::infinity()}},
            fixed, free);
    });
    require_invalid([&] {
        run(c, {{10, std::numeric_limits<double>::quiet_NaN()}}, free);
    });
    require_invalid([&] {
        run(c, fixed, {{20, std::numeric_limits<double>::infinity()}});
    });
    require_invalid([&] { run(c, {{999, 0.0}}, free); });
    require_invalid([&] { run(c, {{10, 0.0}, {10, 0.0}}, free); });
    require_invalid([&] { run(c, fixed, {{20, 1.0}, {20, 1.0}}); });
    require_invalid([&] { run(c, fixed, {{10, 0.0}}); });
    require_invalid([&] { run(c, fixed, {}); });

    // Isolated two-electrode Maxwell matrix has a constant-voltage nullspace.
    require_invalid([&] {
        run(c, {}, {{10, -20.0}, {20, 20.0}});
    });
    // A badly conditioned floating block is not silently inverted.
    require_invalid([&] {
        run({{10, 20}, {1.0, 0.0, 0.0, 1.0e-14}}, {},
            {{10, 1.0}, {20, 1.0}});
    });
}
