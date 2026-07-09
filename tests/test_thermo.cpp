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
    EXPECT_DOUBLE_EQ(state.xk10, 4.199e-17);
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
