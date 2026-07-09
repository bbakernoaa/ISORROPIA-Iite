#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(ThermodynamicsTest, EquilibriumConstantsAtStandardTemp) {
    Isorropia::State state;
    state.temp = 298.15; // standard K
    state.calculate_equilibrium_constants();
    
    // Assert known values at standard temperature (298.15 K)
    EXPECT_DOUBLE_EQ(state.xk1, 1.015e-2);
    EXPECT_DOUBLE_EQ(state.xk21, 57.639);
    EXPECT_DOUBLE_EQ(state.xk22, 1.805e-5);
    EXPECT_DOUBLE_EQ(state.xk3, 1.971e6);
    EXPECT_DOUBLE_EQ(state.xk4, 2.511e6);
    EXPECT_DOUBLE_EQ(state.xk5, 0.4799);
    EXPECT_DOUBLE_EQ(state.xk7, 1.817);
    EXPECT_DOUBLE_EQ(state.xk8, 37.661);
    EXPECT_DOUBLE_EQ(state.xk10, 5.746e-17);
    EXPECT_DOUBLE_EQ(state.xkw, 1.010e-14);

    // Assert derived variables are correct
    EXPECT_DOUBLE_EQ(state.xk2, state.xk21 * state.xk22);
    EXPECT_DOUBLE_EQ(state.xk42, state.xk4 / state.xk41);
    EXPECT_DOUBLE_EQ(state.xk32, state.xk3 / state.xk31);
}

TEST(ThermodynamicsTest, EquilibriumConstantsTemperatureScaling) {
    Isorropia::State state_std;
    state_std.temp = 298.15;
    state_std.calculate_equilibrium_constants();

    Isorropia::State state_cold;
    state_cold.temp = 280.15;
    state_cold.calculate_equilibrium_constants();

    // Check that temperature dependence scales correctly
    // e.g. xk1 should decrease or change according to van 't Hoff scaling
    EXPECT_NE(state_std.xk1, state_cold.xk1);
    EXPECT_NE(state_std.xk21, state_cold.xk21);
    EXPECT_NE(state_std.xkw, state_cold.xkw);
}

TEST(ThermodynamicsTest, ZSRWaterUptake) {
    Isorropia::State state;
    state.initialize_constants();
    state.initialize_water_activities();
    
    // Set meteorological conditions and Case info
    state.rh = 0.80; // 80% RH
    state.scase = "A"; // Case A Metastable Liquid NH4-SO4
    
    // Set liquid concentrations: SO4-- and NH4+ (which CALCMR will map to active salts)
    state.molal[5] = 1.014e-2; // SO4--
    state.molal[6] = 6.487e-5; // HSO4-
    
    // Set organic concentrations: Org = 10 ug/m3, k_org = 0.15, rho_org = 1.0 g/cm3 (1000 kg/m3)
    state.org[0] = 10.0;
    state.org[1] = 0.15;
    state.org[2] = 1000.0;

    // Call ZSR calculator
    state.cal_cmr();

    // Verify water components and total water
    EXPECT_GT(state.water, 0.0);
    EXPECT_GT(state.watcmp[23], 0.0); // organic water contribution must be positive
    EXPECT_NEAR(state.watcmp[23], 6.0, 1e-3); // WatOrg should be exactly 6.0 under these parameters
    EXPECT_GT(state.watcmp[3], 0.0); // (NH4)2SO4 water uptake
}
