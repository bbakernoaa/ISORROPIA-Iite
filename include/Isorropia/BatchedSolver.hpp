#pragma once
#include "Isorropia/Solver.hpp"
#include "Isorropia/mdspan.hpp"

namespace Isorropia {

inline void solve_batched_cells(double* raw_inputs, double* raw_outputs, size_t num_cells, double temp, double rh) {
    mdspan<double, 8> inputs(raw_inputs, num_cells);
    mdspan<double, 8> outputs(raw_outputs, num_cells);

    for (size_t i = 0; i < num_cells; ++i) {
        Input cell_input;
        State cell_state;
        
        cell_input.temp = temp;
        cell_input.rh = rh;
        cell_input.iprob = 0; // Forward

        for (size_t k = 0; k < 8; ++k) {
            cell_input.w[k] = inputs(i, k);
        }

        Solver solver;
        solver.solve(cell_input, cell_state);

        for (size_t k = 0; k < 8; ++k) {
            outputs(i, k) = cell_state.waer[k];
        }
    }
}

} // namespace Isorropia
