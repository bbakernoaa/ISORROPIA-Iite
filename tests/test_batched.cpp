#include <gtest/gtest.h>
#include "Isorropia/BatchedSolver.hpp"
#include <vector>

TEST(MDSpanTest, BasicOperations) {
    double data[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    Isorropia::mdspan<double, 3> span(data, 2);

    EXPECT_EQ(span.num_rows(), 2);
    EXPECT_EQ(span.num_cols(), 3);
    EXPECT_DOUBLE_EQ(span(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(span(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(span(0, 2), 3.0);
    EXPECT_DOUBLE_EQ(span(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(span(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(span(1, 2), 6.0);

    span(1, 1) = 9.9;
    EXPECT_DOUBLE_EQ(data[4], 9.9);
}

TEST(BatchedSolverTest, SolveThreeCellsParity) {
    const size_t num_cells = 3;
    const double temp = 298.15;
    const double rh = 0.80;

    // Define 3 inputs with 8 components each
    std::vector<double> raw_inputs(num_cells * 8, 0.0);

    // Setup input values
    double sulfate = 1e-6;

    // Cell 0: Sulfate-Poor Case A2
    raw_inputs[0 * 8 + 1] = sulfate;       // w[1] = Sulfate
    raw_inputs[0 * 8 + 2] = 2.10 * sulfate; // w[2] = Ammonia

    // Cell 1: Transition Zone A2_B4_Smooth
    raw_inputs[1 * 8 + 1] = sulfate;
    raw_inputs[1 * 8 + 2] = 2.00 * sulfate;

    // Cell 2: Sulfate-Rich B4
    raw_inputs[2 * 8 + 1] = sulfate;
    raw_inputs[2 * 8 + 2] = 1.50 * sulfate;

    // Allocate output array
    std::vector<double> raw_outputs(num_cells * 8, 0.0);

    // Call batched solver
    Isorropia::solve_batched_cells(raw_inputs.data(), raw_outputs.data(), num_cells, temp, rh);

    // Now run individual solvers to assert exact parity
    for (size_t i = 0; i < num_cells; ++i) {
        Isorropia::Input single_input;
        Isorropia::State single_state;

        single_input.temp = temp;
        single_input.rh = rh;
        single_input.iprob = 0;

        for (size_t k = 0; k < 8; ++k) {
            single_input.w[k] = raw_inputs[i * 8 + k];
        }

        Isorropia::Solver solver;
        solver.solve(single_input, single_state);

        for (size_t k = 0; k < 8; ++k) {
            EXPECT_DOUBLE_EQ(raw_outputs[i * 8 + k], single_state.waer[k]);
        }
    }
}
