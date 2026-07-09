#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(SolverTest, SolveNH4SO4NO3Metastable) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // NH4-SO4-NO3 Metastable Case
    input.w[1] = 1.0; // H2SO4
    input.w[2] = 2.0; // NH3
    input.w[3] = 1.0; // HNO3
    input.org[0] = 10.0; // Org
    input.rh = 0.80; // 80% RH
    input.temp = 298.15; // standard Temp

    solver.solve(input, state);

    // Assert that constants and properties are loaded into state
    EXPECT_DOUBLE_EQ(state.temp, 298.15);
    EXPECT_DOUBLE_EQ(state.rh, 0.80);
    EXPECT_DOUBLE_EQ(state.imw[0], 23.0);
    EXPECT_DOUBLE_EQ(state.drnh42s4, 0.7997);

    // Verify solver outputs mapped to standard speciation defaults (Run 1)
    EXPECT_NEAR(state.water, 8.088, 1e-4);
    EXPECT_NEAR(state.gnh3, 1.595, 1e-4);
    EXPECT_NEAR(state.ghno3, 0.7805, 1e-4);
    EXPECT_NEAR(state.molal[1], 2.170e-5, 1e-8); // H+
}
