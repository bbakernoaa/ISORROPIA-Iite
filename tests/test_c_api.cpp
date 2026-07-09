#include <gtest/gtest.h>
#include "Isorropia/Isorropia.h"
#include <cmath>

TEST(CAPITest, SolveCCompatible) {
    IsorropiaInput input = {0.0};
    IsorropiaState state = {0.0};

    // Set up standard NH4-SO4-NO3 metastable case (Run 1)
    input.w[1] = 1.0;  // SO4
    input.w[2] = 2.0;  // NH3
    input.w[3] = 1.0;  // HNO3
    input.org[0] = 10.0; // Org
    input.org[1] = 0.15; // k_org
    input.org[2] = 1000.0; // density
    input.rh = 0.80;   // 80% RH
    input.temp = 298.15; // 298.15K
    input.iprob = 0;   // Forward Problem

    // Call standard C solve entry point
    isorropia_solve_c(&input, &state);

    // Verify linkage values and calculations matching C++ expectations
    EXPECT_DOUBLE_EQ(state.temp, 298.15);
    EXPECT_DOUBLE_EQ(state.rh, 0.80);
    EXPECT_EQ(state.num_errors, 0);

    // Verify outputs have exact match
    EXPECT_NEAR(state.water, 8.088, 1e-4);
    EXPECT_NEAR(state.gnh3, 1.595, 1e-4);
    EXPECT_NEAR(state.ghno3, 0.7805, 1e-4);
    EXPECT_NEAR(state.molal[1], 2.170e-5, 1e-8); // H+
}
