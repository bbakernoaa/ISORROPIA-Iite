#include "Isorropia/Solver.hpp"

namespace Isorropia {

Solver::Solver() = default;

void Solver::solve(const Input& input, State& state) {
    state.clear_errors();

    // 1. Replicate INIT: Copy input variables to State COMMON equivalents
    state.temp = input.temp;
    state.rh   = input.rh;
    state.w    = input.w;
    state.waer = input.waer;
    state.org  = input.org;

    // 2. Initialize constants, lookup tables, deliquescence coefficients
    state.initialize_constants();
    state.initialize_water_activities();
    state.initialize_drh();
    state.calculate_equilibrium_constants();

    // 3. Dispatch based on formulation type (0 = Forward, 1 = Reverse)
    double nitrate_sum = input.w[3] + input.waer[3]; // HNO3 or aerosol nitrate
    double crustal_sum = input.w[5] + input.w[6] + input.w[7] + input.waer[5] + input.waer[6] + input.waer[7]; // Ca, K, Mg
    double marine_sum  = input.w[0] + input.w[4] + input.waer[0] + input.waer[4]; // Na, Cl
    
    if (input.iprob == 1) {
        // Reverse problem routing
        if (crustal_sum > state.tiny) {
            state.push_error(104, "Case 4 Reverse is not implemented yet in the port skeleton");
        } else if (marine_sum > state.tiny) {
            state.push_error(103, "Case 3 Reverse is not implemented yet in the port skeleton");
        } else if (nitrate_sum > state.tiny) {
            isrp2r(input, state);
        } else {
            isrp1r(input, state);
        }
    } else {
        // Forward problem routing
        if (crustal_sum > state.tiny) {
            isrp4f(input, state);
        } else if (marine_sum > state.tiny) {
            isrp3f(input, state);
        } else if (nitrate_sum > state.tiny) {
            isrp2f(input, state);
        } else {
            isrp1f(input, state);
        }
    }
}

} // namespace Isorropia
