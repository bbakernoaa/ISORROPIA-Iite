#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(SolverTest, SolveNH4SO4NO3Metastable) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // NH4-SO4-NO3 Metastable Case (Run 1 inputs scaled to standard internal mol/m3, kg/m3)
    input.w[1] = (1.0 / 98.0) * 1e-6; // Sulfate: 1.020e-8 mol/m3
    input.w[2] = (2.0 / 17.0) * 1e-6; // Ammonia: 1.176e-7 mol/m3
    input.w[3] = (1.0 / 63.0) * 1e-6; // Nitrate: 1.587e-8 mol/m3
    
    input.org[0] = 10.0 * 1e-9;        // Organic mass: 10.0e-9 kg/m3
    input.org[1] = 0.15;               // Organic kappa
    input.org[2] = 1.0 * 1000.0;       // Organic density: 1000.0 kg/m3

    input.rh = 0.80; // 80% RH
    input.temp = 298.15; // standard Temp

    solver.solve(input, state);

    // Assert that constants and properties are loaded into state
    EXPECT_DOUBLE_EQ(state.temp, 298.15);
    EXPECT_DOUBLE_EQ(state.rh, 0.80);
    EXPECT_DOUBLE_EQ(state.imw[0], 23.0);
    EXPECT_DOUBLE_EQ(state.drnh42s4, 0.7997);

    // Verify solver outputs mapped to standard physical scales
    EXPECT_NEAR(state.water, 7.89625e-9, 1e-12);
    EXPECT_NEAR(state.gnh3 * 1e6, 9.5772e-2, 1e-4);  // Gaseous Ammonia in umol/m3
    EXPECT_NEAR(state.ghno3 * 1e6, 1.4363e-2, 1e-4); // Gaseous Nitric Acid in umol/m3
    EXPECT_NEAR(state.molal[1], 1.994e-12, 1e-14);  // H+ molality (mol/kg)
}
