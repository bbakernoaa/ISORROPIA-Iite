#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(CrustalTest, ForwardCrustalSolve) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;
    
    // Set up a standard crustal configuration (Case 4)
    input.w[0] = 1.0; // Na
    input.w[1] = 2.0; // H2SO4
    input.w[2] = 2.0; // NH3
    input.w[3] = 1.0; // HNO3
    input.w[5] = 0.5; // Ca
    input.w[6] = 0.2; // K
    input.w[7] = 0.3; // Mg
    input.rh = 0.80; // 80% RH
    input.temp = 298.15; // standard K

    solver.solve(input, state);

    // Verify appropriate routing and initial speciation setup
    EXPECT_EQ(state.scase, "4F"); // Route to Case 4 Forward Crustal
    EXPECT_GT(state.water, 0.0); // ZSR must resolve water content
    EXPECT_EQ(state.num_errors, 0); // No error diagnostics triggered
}
