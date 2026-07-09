#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(ReverseTest, SolveNH4SO4Reverse) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // Set up standard NH4-SO4 system in Reverse mode (iprob = 1)
    input.waer[1] = 1.0; // Aerosol H2SO4
    input.waer[2] = 2.0; // Aerosol NH3
    input.iprob = 1; // Reverse Problem!
    input.rh = 0.80;

    solver.solve(input, state);

    // Verify appropriate routing and case setup
    EXPECT_EQ(state.scase, "S2"); // Route to Case 1 Reverse (routes specifically to S2 subcase)
    EXPECT_EQ(state.num_errors, 0);
}

TEST(ReverseTest, SolveNH4SO4NO3Reverse) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // Set up NH4-SO4-NO3 system in Reverse mode (iprob = 1)
    input.waer[1] = 1.0; // Aerosol H2SO4
    input.waer[2] = 2.0; // Aerosol NH3
    input.waer[3] = 1.0; // Aerosol HNO3
    input.iprob = 1; // Reverse Problem!
    input.rh = 0.80;

    solver.solve(input, state);

    // Verify appropriate routing and case setup
    EXPECT_EQ(state.scase, "N3"); // Route to Case 2 Reverse (routes specifically to N3 subcase)
    EXPECT_EQ(state.num_errors, 0);
}

TEST(ReverseTest, SolveMarineReverse) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // Set up Case 3 Marine reverse
    input.waer[0] = 1.0; // Aerosol Na
    input.waer[1] = 1.0; // Aerosol H2SO4
    input.waer[4] = 1.0; // Aerosol Cl
    input.iprob = 1; // Reverse!
    input.rh = 0.80;

    solver.solve(input, state);

    EXPECT_EQ(state.scase, "I6"); // Route to Case 3 Reverse (routes to I6 subcase)
    EXPECT_GT(state.water, 0.0);
    EXPECT_EQ(state.num_errors, 0);
}

TEST(ReverseTest, SolveCrustalReverse) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;

    // Set up Case 4 Crustal reverse
    input.waer[0] = 1.0; // Aerosol Na
    input.waer[1] = 2.0; // Aerosol H2SO4
    input.waer[5] = 0.5; // Aerosol Ca
    input.iprob = 1; // Reverse!
    input.rh = 0.80;

    solver.solve(input, state);

    EXPECT_EQ(state.scase, "4R"); // Route to Case 4 Reverse
    EXPECT_GT(state.water, 0.0);
    EXPECT_EQ(state.num_errors, 0);
}
