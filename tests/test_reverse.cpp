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

TEST(ReverseTest, SolveMarineReverse) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // Set up Case 3 Marine reverse
    input.w[0] = 1.0; // Na
    input.w[1] = 1.0; // H2SO4
    input.w[4] = 1.0; // Cl
    input.iprob = 1; // Reverse!
    input.rh = 0.80;

    solver.solve(input, state);

    EXPECT_EQ(state.scase, "3R"); // Route to Case 3 Reverse
    EXPECT_GT(state.water, 0.0);
    EXPECT_EQ(state.num_errors, 0);
}

TEST(ReverseTest, SolveCrustalReverse) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // Set up Case 4 Crustal reverse
    input.w[0] = 1.0; // Na
    input.w[1] = 2.0; // H2SO4
    input.w[5] = 0.5; // Ca
    input.iprob = 1; // Reverse!
    input.rh = 0.80;

    solver.solve(input, state);

    EXPECT_EQ(state.scase, "4R"); // Route to Case 4 Reverse
    EXPECT_GT(state.water, 0.0);
    EXPECT_EQ(state.num_errors, 0);
}
