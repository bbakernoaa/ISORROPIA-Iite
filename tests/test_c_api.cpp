#include <gtest/gtest.h>
#include "Isorropia/Isorropia.h"
#include <cmath>

TEST(CAPITest, SolveCCompatible) {
    IsorropiaInput input = {0.0};
    IsorropiaState state = {0.0};

    // Set up standard NH4-SO4-NO3 metastable case (Run 1 inputs scaled to standard mol/m3)
    input.w[1] = (1.0 / 98.0) * 1e-6;  // H2SO4 -> mol/m3
    input.w[2] = (2.0 / 17.0) * 1e-6;  // NH3 -> mol/m3
    input.w[3] = (1.0 / 63.0) * 1e-6;  // HNO3 -> mol/m3
    input.org[0] = 10.0 * 1e-9;        // Org -> kg/m3
    input.org[1] = 0.15;               // k_org
    input.org[2] = 1.0 * 1000.0;       // density -> kg/m3
    input.rh = 0.80;   // 80% RH
    input.temp = 298.15; // 298.15K
    input.iprob = 0;   // Forward Problem

    // Call standard C solve entry point
    isorropia_solve_c(&input, &state);

    // Verify linkage values and calculations matching C++ expectations
    EXPECT_DOUBLE_EQ(state.temp, 298.15);
    EXPECT_DOUBLE_EQ(state.rh, 0.80);
    EXPECT_EQ(state.num_errors, 0);

    // Verify outputs have exact match in standard physical units
    EXPECT_NEAR(state.water, 8.083366e-9, 1e-11);
    EXPECT_NEAR(state.gnh3 * 1e6, 9.3891e-2, 1e-4);  // Gaseous Ammonia in umol/m3
    EXPECT_NEAR(state.ghno3 * 1e6, 1.2439e-2, 1e-4); // Gaseous Nitric Acid in umol/m3
    EXPECT_NEAR(state.molal[1], 2.0889e-11, 1e-13);  // H+ molality (mol/kg)
}
