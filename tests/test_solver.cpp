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
    EXPECT_NEAR(state.water, 8.083366e-9, 1e-11);
    EXPECT_NEAR(state.gnh3 * 1e6, 9.3891e-2, 1e-4);  // Gaseous Ammonia in umol/m3
    EXPECT_NEAR(state.ghno3 * 1e6, 1.2439e-2, 1e-4); // Gaseous Nitric Acid in umol/m3
    EXPECT_NEAR(state.molal[1], 2.170486e-11, 1e-13);  // H+ molality (mol/kg)
}

TEST(SolverTest, Case1SmoothTransitions) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // We want only Sulfate and Ammonia present, no Nitrate, no Sodium, no Crustals, no Chlorine.
    // Set standard temperature and relative humidity
    input.temp = 298.15;
    input.rh = 0.80;
    input.iprob = 0; // Forward solver

    // Let's test the five regions by varying Ammonia (input.w[2]) for a fixed Sulfate (input.w[1] = 1e-6)
    double sulfate = 1e-6;
    input.w[1] = sulfate;

    // Case A: Pure Sulfate-Poor (A2)
    // sulrat > 2.05, e.g. 2.10 -> Ammonia = 2.10e-6
    input.w[2] = 2.10 * sulfate;
    solver.solve(input, state);
    EXPECT_EQ(state.scase, "A2");

    // Case B: A2-B4 Transition Zone (A2_B4_Smooth)
    // sulrat in [1.95, 2.05], e.g. 2.00 -> Ammonia = 2.00e-6
    input.w[2] = 2.00 * sulfate;
    solver.solve(input, state);
    EXPECT_EQ(state.scase, "A2_B4_Smooth");
    EXPECT_GT(state.water, 0.0);

    // Case C: Pure Sulfate-Rich / No Free Acid (B4)
    // sulrat in (1.05, 1.95), e.g. 1.50 -> Ammonia = 1.50e-6
    input.w[2] = 1.50 * sulfate;
    solver.solve(input, state);
    EXPECT_EQ(state.scase, "B4");

    // Case D: B4-C2 Transition Zone (B4_C2_Smooth)
    // sulrat in [0.95, 1.05], e.g. 1.00 -> Ammonia = 1.00e-6
    input.w[2] = 1.00 * sulfate;
    solver.solve(input, state);
    EXPECT_EQ(state.scase, "B4_C2_Smooth");

    // Case E: Pure Sulfate-Rich / Free Acid (C2)
    // sulrat < 0.95, e.g. 0.80 -> Ammonia = 0.80e-6
    input.w[2] = 0.80 * sulfate;
    solver.solve(input, state);
    EXPECT_EQ(state.scase, "C2");
}

TEST(SolverTest, Case2SmoothTransitions) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    input.temp = 298.15;
    input.rh = 0.80;
    input.iprob = 0; // Forward solver

    double sulfate = 1e-6;
    input.w[1] = sulfate;
    input.w[3] = 0.5e-6; // Nitrate present

    // Case A: Pure Sulfate-Poor (D3)
    input.w[2] = 2.10 * sulfate;
    solver.solve(input, state);
    EXPECT_EQ(state.scase, "D3");

    // Case B: D3-E4 Transition Zone (D3_E4_Smooth)
    input.w[2] = 2.00 * sulfate;
    solver.solve(input, state);
    EXPECT_EQ(state.scase, "D3_E4_Smooth");
    EXPECT_GT(state.water, 0.0);

    // Case C: Pure Sulfate-Rich / No Free Acid (E4)
    input.w[2] = 1.50 * sulfate;
    solver.solve(input, state);
    EXPECT_EQ(state.scase, "E4");

    // Case D: E4-F2 Transition Zone (E4_F2_Smooth)
    input.w[2] = 1.00 * sulfate;
    solver.solve(input, state);
    EXPECT_EQ(state.scase, "E4_F2_Smooth");

    // Case E: Pure Sulfate-Rich / Free Acid (F2)
    input.w[2] = 0.80 * sulfate;
    solver.solve(input, state);
    EXPECT_EQ(state.scase, "F2");
}
