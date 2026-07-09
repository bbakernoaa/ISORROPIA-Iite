#include "Isorropia/Solver.hpp"

namespace Isorropia {

Solver::Solver() = default;

void Solver::solve(const Input& input, State& state) {
    state.clear_errors();

    // 1. Replicate INIT: Copy input meteorological variables to State COMMON equivalents
    state.temp = input.temp;
    state.rh   = input.rh;

    // 2. Initialize constants, lookup tables, deliquescence coefficients
    state.initialize_constants();
    state.initialize_water_activities();
    state.initialize_drh();
    state.calculate_equilibrium_constants();

    // 3. Dispatch to specific forward solver case.
    // If Nitrate/HNO3 are non-zero, it is Case 2 (NH4-SO4-NO3 Metastable system).
    // Otherwise, it is Case 1 (NH4-SO4-H2O Metastable system).
    double nitrate_sum = input.w[3] + input.waer[3]; // HNO3 or aerosol nitrate
    
    if (nitrate_sum > state.tiny) {
        isrp2f(input, state);
    } else {
        isrp1f(input, state);
    }
}

} // namespace Isorropia
