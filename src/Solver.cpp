#include "Isorropia/Solver.hpp"
#include "Isorropia/Isorropia.h"
#include <algorithm>

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

    // 2b. Replicate Fortran INIT1-4 flooring of total/aerosol concentrations to
    // TINY (absent species become 1e-20, NOT 0). This is required for bit-faithful
    // regime-boundary selection (e.g. SULRATW vs SO4RAT ties). Forward floors W;
    // reverse floors WAER and zeros W. The dispatch sums above use `input` (the
    // original, unfloored values), so routing is unaffected.
    if (input.iprob == 1) {
        for (auto& v : state.waer) v = std::max(v, state.tiny);
        state.w.fill(0.0);
    } else {
        for (auto& v : state.w) v = std::max(v, state.tiny);
    }

    // 3. Dispatch based on formulation type (0 = Forward, 1 = Reverse)
    double nitrate_sum = input.w[3] + input.waer[3]; // HNO3 or aerosol nitrate
    double crustal_sum = input.w[5] + input.w[6] + input.w[7] + input.waer[5] + input.waer[6] + input.waer[7]; // Ca, K, Mg
    double marine_sum  = input.w[0] + input.w[4] + input.waer[0] + input.waer[4]; // Na, Cl

    if (input.iprob == 1) {
        // Reverse problem routing
        if (crustal_sum > state.tiny) {
            isrp4r(input, state);
        } else if (marine_sum > state.tiny) {
            isrp3r(input, state);
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

// Standard C Linkage Implementation
extern "C" {

void isorropia_solve_c(const IsorropiaInput* input, IsorropiaState* state) {
    if (!input || !state) return;

    // 1. Map standard flat C Input struct into thread-local C++ Input structure
    Isorropia::Input cpp_input;
    std::copy(std::begin(input->w), std::end(input->w), cpp_input.w.begin());
    std::copy(std::begin(input->org), std::end(input->org), cpp_input.org.begin());
    std::copy(std::begin(input->waer), std::end(input->waer), cpp_input.waer.begin());
    cpp_input.temp  = input->temp;
    cpp_input.rh    = input->rh;
    cpp_input.iprob = input->iprob;
    cpp_input.nadj  = input->nadj;

    // 2. Instantiate local C++ State
    Isorropia::State cpp_state;

    // 3. Call the high-performance C++ solver
    Isorropia::Solver solver;
    solver.solve(cpp_input, cpp_state);

    // 4. Map C++ State output fields back into flat C State structure
    state->temp = cpp_state.temp;
    state->rh   = cpp_state.rh;
    std::copy(cpp_state.w.begin(), cpp_state.w.end(), state->w);
    std::copy(cpp_state.waer.begin(), cpp_state.waer.end(), state->waer);
    std::copy(cpp_state.org.begin(), cpp_state.org.end(), state->org);

    std::copy(cpp_state.molal.begin(), cpp_state.molal.end(), state->molal);
    for (size_t i = 0; i < 23; ++i) {
        state->molalr[i] = cpp_state.molalr[i];
        state->gama[i]   = cpp_state.gama[i];
        state->zz[i]     = cpp_state.zz[i];
        state->gamou[i]  = cpp_state.gamou[i];
        state->gamin[i]  = cpp_state.gamin[i];
        state->m0[i]     = cpp_state.m0[i];
    }
    std::copy(cpp_state.z.begin(), cpp_state.z.end(), state->z);
    std::copy(cpp_state.gasaq.begin(), cpp_state.gasaq.end(), state->gasaq);

    state->actmod = cpp_state.actmod;
    state->epsact = cpp_state.epsact;
    state->coh    = cpp_state.coh;
    state->chno3  = cpp_state.chno3;
    state->chcl   = cpp_state.chcl;
    state->water  = cpp_state.water;
    state->ionic  = cpp_state.ionic;
    std::copy(cpp_state.watcmp.begin(), cpp_state.watcmp.end(), state->watcmp);

    state->frst   = cpp_state.frst ? 1 : 0;
    state->calain = cpp_state.calain ? 1 : 0;
    state->calaou = cpp_state.calaou ? 1 : 0;
    state->dryf   = cpp_state.dryf ? 1 : 0;

    // Map solid salt concentrations
    state->ch2so4  = cpp_state.ch2so4;
    state->cnh42s4 = cpp_state.cnh42s4;
    state->cnh4hs4 = cpp_state.cnh4hs4;
    state->cnacl   = cpp_state.cnacl;
    state->cna2so4 = cpp_state.cna2so4;
    state->cnano3  = cpp_state.cnano3;
    state->cnh4no3 = cpp_state.cnh4no3;
    state->cnh4cl  = cpp_state.cnh4cl;
    state->cnahso4 = cpp_state.cnahso4;
    state->clc     = cpp_state.clc;
    state->ccaso4  = cpp_state.ccaso4;
    state->ccano32 = cpp_state.ccano32;
    state->ccacl2  = cpp_state.ccacl2;
    state->ck2so4  = cpp_state.ck2so4;
    state->ckhso4  = cpp_state.ckhso4;
    state->ckno3   = cpp_state.ckno3;
    state->ckcl    = cpp_state.ckcl;
    state->cmgso4  = cpp_state.cmgso4;
    state->cmgno32 = cpp_state.cmgno32;
    state->cmgcl2  = cpp_state.cmgcl2;

    // Map gas concentrations
    state->gnh3  = cpp_state.gnh3;
    state->ghno3 = cpp_state.ghno3;
    state->ghcl  = cpp_state.ghcl;

    // Diagnostic counts
    state->num_errors = static_cast<int>(cpp_state.num_errors);
}

}
