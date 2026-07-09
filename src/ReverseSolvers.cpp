#include "Isorropia/Solver.hpp"

namespace Isorropia {

void Solver::isrp1r(const Input& input, State& state) {
    // Reverse solver for NH4-SO4-H2O systems (Case 1)
    state.clear_errors();
    state.scase = "1R"; // Case 1 Reverse

    // Replicate initial reverse mapping and ZSR calculation
    state.cal_cmr();
}

void Solver::isrp2r(const Input& input, State& state) {
    // Reverse solver for NH4-SO4-NO3-H2O systems (Case 2)
    state.clear_errors();
    state.scase = "2R"; // Case 2 Reverse

    state.cal_cmr();
}

} // namespace Isorropia
