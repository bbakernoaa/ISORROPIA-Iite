#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(StateTest, ErrorStackPushAndClear) {
    Isorropia::State state;
    
    state.push_error(1, "Test Error 1");
    state.push_error(2, "Test Error 2");
    
    EXPECT_EQ(state.num_errors, 2);
    EXPECT_FALSE(state.stack_overflow);
    EXPECT_EQ(state.error_stack[0].code, 1);
    EXPECT_EQ(state.error_stack[0].message, "Test Error 1");
    EXPECT_EQ(state.error_stack[1].code, 2);
    
    state.clear_errors();
    EXPECT_EQ(state.num_errors, 0);
    EXPECT_FALSE(state.stack_overflow);
}

TEST(StateTest, ErrorStackOverflow) {
    Isorropia::State state;
    
    for (int i = 0; i < 26; ++i) {
        state.push_error(i, "Overflow Test");
    }
    
    EXPECT_EQ(state.num_errors, 25);
    EXPECT_TRUE(state.stack_overflow);
}

TEST(StateTest, InitializeConstantsAndWeights) {
    Isorropia::State state;
    state.initialize_constants();
    
    EXPECT_DOUBLE_EQ(state.r, 82.0567e-6);
    EXPECT_DOUBLE_EQ(state.imw[0], 23.0); // Na+
    EXPECT_DOUBLE_EQ(state.imw[2], 18.0); // NH4+
    EXPECT_DOUBLE_EQ(state.wmw[1], 98.0); // H2SO4
    EXPECT_DOUBLE_EQ(state.smw[5], 132.0); // (NH4)2SO4
}

TEST(StateTest, WaterActivityGrids) {
    Isorropia::State state;
    state.initialize_water_activities();
    
    // Assert copy fidelity for known parsed tables (e.g. awab and awac)
    EXPECT_NEAR(state.awab[0], 3.128400e+02, 1e-6);
    EXPECT_NEAR(state.awab[99], 1.000000e-01, 1e-6);
    EXPECT_NEAR(state.awac[0], 1.209000e+03, 1e-6);
    EXPECT_NEAR(state.awac[99], 1.000000e-01, 1e-6);
    EXPECT_NEAR(state.awcs[0], 0.0, 1e-15); // should be filled with zero
}
