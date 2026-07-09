#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(ReverseTest, SolveNH4SO4Reverse) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // Set up standard NH4-SO4 system in Reverse mode (iprob = 1)
    input.w[1] = 1.0; // H2SO4
    input.w[2] = 2.0; // NH3
    input.iprob = 1; // Reverse Problem!
    input.rh = 0.80;

    solver.solve(input, state);

    // Verify appropriate routing and case setup
    EXPECT_EQ(state.scase, "1R"); // Route to Case 1 Reverse
    EXPECT_EQ(state.num_errors, 0);
}

TEST(ReverseTest, SolveNH4SO4NO3Reverse) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // Set up NH4-SO4-NO3 system in Reverse mode (iprob = 1)
    input.w[1] = 1.0; // H2SO4
    input.w[2] = 2.0; // NH3
    input.w[3] = 1.0; // HNO3
    input.iprob = 1; // Reverse Problem!
    input.rh = 0.80;

    solver.solve(input, state);

    // Verify appropriate routing and case setup
    EXPECT_EQ(state.scase, "2R"); // Route to Case 2 Reverse
    EXPECT_EQ(state.num_errors, 0);
}
