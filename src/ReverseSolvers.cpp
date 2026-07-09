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

void Solver::isrp3r(const Input& input, State& state) {
    // Reverse solver for Na-NH4-SO4-NO3-Cl-H2O marine systems (Case 3)
    state.clear_errors();
    state.scase = "3R"; // Case 3 Reverse

    state.cal_cmr();
}

void Solver::isrp4r(const Input& input, State& state) {
    // Reverse solver for Na-NH4-SO4-NO3-Cl-Ca-K-Mg-H2O crustal systems (Case 4)
    state.clear_errors();
    state.scase = "4R"; // Case 4 Reverse

    state.cal_cmr();
}

} // namespace Isorropia
