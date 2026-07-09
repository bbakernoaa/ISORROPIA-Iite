#include "Isorropia/Solver.hpp"
#include <algorithm>
#include <cmath>

namespace Isorropia {

//=======================================================================
// REVERSE SOLVER GAS PHASE RECONSTRUCTIONS (CALCNH3P & CALCNHP equivalents)
//=======================================================================
void Solver::cal_cnh3p(State& state) {
    if (state.water <= state.tiny || state.molal[1] <= state.tiny) {
        state.gnh3 = state.tiny;
        return;
    }

    // A1 = (XK2/XKW) * R * TEMP * (GAMA(10)/GAMA(5))^2
    double a1 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
    state.gnh3 = std::max(state.tiny, state.molal[2] / state.molal[1] / a1); // NH4+ / H+ / A1
}

void Solver::cal_cnhp(const Input& input, State& state) {
    if (state.water <= state.tiny) {
        state.ghcl  = std::max(state.w[4] - state.molal[4], state.tiny);
        state.ghno3 = std::max(state.w[3] - state.molal[3], state.tiny);
        return;
    }

    // A3 = XK3 * R * TEMP * (WATER / GAMA(11))^2 -> HCl
    double a3 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);
    // A4 = XK4 * R * TEMP * (WATER / GAMA(10))^2 -> HNO3
    double a4 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);

    // Temp increase H+ due to added volatile acids
    state.molal[1] += (state.waer[3] + state.waer[4]);

    double delt = 0.0;
    // Assume del_t HNO3 >> del_t HCl
    cal_niaq(state.waer[3], state.molal[1] + state.molal[3] + state.molal[4], delt, state);
    state.molal[1] -= delt;
    state.molal[3]  = std::max(state.tiny, state.waer[3] - delt); // NO3-
    state.gasaq[2]  = delt; // HNO3aq -> index 2

    cal_claq(state.waer[4], state.molal[1] + state.molal[3] + state.molal[4], delt, state);
    state.molal[1] -= delt;
    state.molal[4]  = std::max(state.tiny, state.waer[4] - delt); // Cl-
    state.gasaq[1]  = delt; // HClaq -> index 1

    state.ghno3 = std::max(state.tiny, state.molal[1] * state.molal[3] / a4);
    state.ghcl  = std::max(state.tiny, state.molal[1] * state.molal[4] / a3);
}

//=======================================================================
// REVERSE PROBLEM SOLVERS
//=======================================================================
void Solver::isrp1r(const Input& input, State& state) {
    state.clear_errors();
    state.scase = "1R";
    state.actmod = 1;

    // 1. Map input aerosol concentrations directly to liquid ion concentrations
    state.molal[0] = state.waer[0]; // Na+ -> maps to molal[0]
    state.molal[2] = state.waer[2]; // NH4+ -> maps to molal[2]
    state.molal[5] = state.waer[1]; // SO4-- -> maps to molal[5]
    state.molal[6] = 0.0;           // HSO4- -> maps to molal[6]

    // 2. Charge balance pH evaluation (MOLAL(1) = H+)
    double smin = 2.0 * state.molal[5] - state.molal[0] - state.molal[2];
    double hi = 0.0, ohi = 0.0;
    cal_cph(smin, hi, ohi, state);
    state.molal[1] = hi;

    // 3. Compute water uptake and activities
    state.cal_cmr();
    state.cal_act1();

    // 4. Calculate Gas partitionings
    cal_cnh3p(state);
}

void Solver::isrp2r(const Input& input, State& state) {
    state.clear_errors();
    state.scase = "2R";
    state.actmod = 2;

    // 1. Map input aerosol concentrations
    state.molal[0] = state.waer[0]; // Na+
    state.molal[2] = state.waer[2]; // NH4+
    state.molal[3] = state.waer[3]; // NO3-
    state.molal[5] = state.waer[1]; // SO4--
    state.molal[6] = 0.0;           // HSO4-

    // 2. Charge balance pH
    double smin = 2.0 * state.molal[5] + state.molal[3] - state.molal[0] - state.molal[2];
    double hi = 0.0, ohi = 0.0;
    cal_cph(smin, hi, ohi, state);
    state.molal[1] = hi;

    // 3. Compute water uptake and activities
    state.cal_cmr();
    state.cal_act2();

    // 4. Gas partitionings
    cal_cnh3p(state);
    cal_cnhp(input, state);
}

void Solver::isrp3r(const Input& input, State& state) {
    state.clear_errors();
    state.scase = "3R";
    state.actmod = 3;

    // 1. Map input aerosol concentrations
    state.molal[0] = state.waer[0]; // Na+
    state.molal[2] = state.waer[2]; // NH4+
    state.molal[3] = state.waer[3]; // NO3-
    state.molal[4] = state.waer[4]; // Cl-
    state.molal[5] = state.waer[1]; // SO4--
    state.molal[6] = 0.0;           // HSO4-

    // 2. Charge balance pH
    double smin = 2.0 * state.molal[5] + state.molal[3] + state.molal[4] - state.molal[0] - state.molal[2];
    double hi = 0.0, ohi = 0.0;
    cal_cph(smin, hi, ohi, state);
    state.molal[1] = hi;

    // 3. Compute water uptake and activities
    state.cal_cmr();
    state.cal_act2(); // Uses Standard Case 2/3 active pairs

    // 4. Gas partitionings
    cal_cnh3p(state);
    cal_cnhp(input, state);
}

void Solver::isrp4r(const Input& input, State& state) {
    state.clear_errors();
    state.scase = "4R";
    state.actmod = 4;

    // 1. Map input aerosol concentrations (replicates INIT4)
    state.molal[0] = state.waer[0]; // Na+
    state.molal[2] = state.waer[2]; // NH4+
    state.molal[3] = state.waer[3]; // NO3-
    state.molal[4] = state.waer[4]; // Cl-
    state.molal[5] = state.waer[1]; // SO4--
    state.molal[6] = 0.0;           // HSO4-
    state.molal[7] = state.waer[5]; // Ca++
    state.molal[8] = state.waer[6]; // K+
    state.molal[9] = state.waer[7]; // Mg++

    // 2. Charge balance pH including crustal charges
    double smin = 2.0 * state.molal[5] + state.molal[3] + state.molal[4] 
                - state.molal[0] - state.molal[2] 
                + 2.0 * state.molal[7] + state.molal[8] + 2.0 * state.molal[9];
    double hi = 0.0, ohi = 0.0;
    cal_cph(smin, hi, ohi, state);
    state.molal[1] = hi;

    // 3. Compute water uptake and activities
    state.cal_cmr();
    state.cal_act2();

    // 4. Gas partitionings
    cal_cnh3p(state);
    cal_cnhp(input, state);
}

} // namespace Isorropia
