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

namespace {

static const std::array<double, 14> ASSO4 = {
    1.0E-9, 2.5E-9, 5.0E-9, 7.5E-9, 1.0E-8,
    2.5E-8, 5.0E-8, 7.5E-8, 1.0E-7, 2.5E-7,
    5.0E-7, 7.5E-7, 1.0E-6, 5.0E-6
};

static const std::array<double, 280> ASRAT = {
    1.020464, 0.9998130, 0.9960167, 0.9984423, 1.004004,
    1.010885,  1.018356,  1.026726,  1.034268, 1.043846,
    1.052933,  1.062230,  1.062213,  1.080050, 1.088350,
    1.096603,  1.104289,  1.111745,  1.094662, 1.121594,
    1.268909,  1.242444,  1.233815,  1.232088, 1.234020,
    1.238068,  1.243455,  1.250636,  1.258734, 1.267543,
    1.276948,  1.286642,  1.293337,  1.305592, 1.314726,
    1.323463,  1.333258,  1.343604,  1.344793, 1.355571,
    1.431463,  1.405204,  1.395791,  1.393190, 1.394403,
    1.398107,  1.403811,  1.411744,  1.420560, 1.429990,
    1.439742,  1.449507,  1.458986,  1.468403, 1.477394,
    1.487373,  1.495385,  1.503854,  1.512281, 1.520394,
    1.514464,  1.489699,  1.480686,  1.478187, 1.479446,
    1.483310,  1.489316,  1.497517,  1.506501, 1.515816,
    1.524724,  1.533950,  1.542758,  1.551730, 1.559587,
    1.568343,  1.575610,  1.583140,  1.590440, 1.596481,
    1.567743,  1.544426,  1.535928,  1.533645, 1.535016,
    1.539003,  1.545124,  1.553283,  1.561886, 1.570530,
    1.579234,  1.587813,  1.595956,  1.603901, 1.611349,
    1.618833,  1.625819,  1.632543,  1.639032, 1.645276,
    1.707390,  1.689553,  1.683198,  1.681810, 1.683490,
    1.687477,  1.693148,  1.700084,  1.706917, 1.713507,
    1.719952,  1.726190,  1.731985,  1.737544, 1.742673,
    1.747756,  1.752431,  1.756890,  1.761141, 1.765190,
    1.785657,  1.771851,  1.767063,  1.766229, 1.767901,
    1.771455,  1.776223,  1.781769,  1.787065, 1.792081,
    1.796922,  1.801561,  1.805832,  1.809896, 1.813622,
    1.817292,  1.820651,  1.823841,  1.826871, 1.829745,
    1.822215,  1.810497,  1.806496,  1.805898, 1.807480,
    1.810684,  1.814860,  1.819613,  1.824093, 1.828306,
    1.832352,  1.836209,  1.839748,  1.843105, 1.846175,
    1.849192,  1.851948,  1.854574,  1.857038, 1.859387,
    1.844588,  1.834208,  1.830701,  1.830233, 1.831727,
    1.834665,  1.838429,  1.842658,  1.846615, 1.850321,
    1.853869,  1.857243,  1.860332,  1.863257, 1.865928,
    1.868550,  1.870942,  1.873208,  1.875355, 1.877389,
    1.899556,  1.892637,  1.890367,  1.890165, 1.891317,
    1.893436,  1.896036,  1.898872,  1.901485, 1.903908,
    1.906212,  1.908391,  1.910375,  1.912248, 1.913952,
    1.915621,  1.917140,  1.918576,  1.919934, 1.921220,
    1.928264,  1.923245,  1.921625,  1.921523, 1.922421,
    1.924016,  1.925931,  1.927991,  1.929875, 1.931614,
    1.933262,  1.934816,  1.936229,  1.937560, 1.938769,
    1.939951,  1.941026,  1.942042,  1.943003, 1.943911,
    1.941205,  1.937060,  1.935734,  1.935666, 1.936430,
    1.937769,  1.939359,  1.941061,  1.942612, 1.944041,
    1.945393,  1.946666,  1.947823,  1.948911, 1.949900,
    1.950866,  1.951744,  1.952574,  1.953358, 1.954099,
    1.948985,  1.945372,  1.944221,  1.944171, 1.944850,
    1.946027,  1.947419,  1.948902,  1.950251, 1.951494,
    1.952668,  1.953773,  1.954776,  1.955719, 1.956576,
    1.957413,  1.958174,  1.958892,  1.959571, 1.960213,
    1.977193,  1.975540,  1.975023,  1.975015, 1.975346,
    1.975903,  1.976547,  1.977225,  1.977838, 1.978401,
    1.978930,  1.979428,  1.979879,  1.980302, 1.980686,
    1.981060,  1.981401,  1.981722,  1.982025, 1.982312
};

} // namespace

double Solver::getasr(double so4i, double rhi) {
    // Replicates FUNCTION GETASR and BLOCK DATA AERSR
    double rat = so4i / 1.e-9;
    double a1_val = std::floor(std::log10(rat));
    int ia1 = static_cast<int>(rat / 2.5 / std::pow(10.0, a1_val));

    int inds = static_cast<int>(4.0 * a1_val) + std::min(ia1, 4);
    inds = std::min(std::max(0, inds), 13) + 1; // 1-based index mapping to INDS

    int indr = static_cast<int>(99.0 - rhi * 100.0) + 1;
    indr = std::min(std::max(1, indr), 20);

    int indsl = inds;
    int indsh = std::min(indsl + 1, 14);

    int iposl = (indsl - 1) * 20 + indr - 1; // 0-based offset
    int iposh = (indsh - 1)* 20 + indr - 1;

    double wf = (so4i - ASSO4[indsl - 1]) / (ASSO4[indsh - 1] - ASSO4[indsl - 1] + 1e-7);
    wf = std::min(std::max(wf, 0.0), 1.0);

    return wf * ASRAT[iposh] + (1.0 - wf) * ASRAT[iposl];
}

//=======================================================================
// REVERSE PROBLEM SOLVERS
//=======================================================================
void Solver::isrp1r(const Input& input, State& state) {
    state.clear_errors();
    state.scase = "1R";
    state.actmod = 1;

    double sulratw = 2.0;
    if (state.rh >= state.drnh42s4) {
        sulratw = getasr(state.waer[1], state.rh);
    }
    double sulrat = state.waer[2] / state.waer[1];

    if (sulratw <= sulrat) {
        state.scase = "S2";
        cal_s2(input, state);
    } else if (sulrat >= 1.0) {
        state.w[1] = state.waer[1];
        state.w[2] = state.waer[2];
        state.scase = "B4";
        cal_cb4(input, state);
        cal_cnh3p(state);
    } else {
        state.w[1] = state.waer[1];
        state.w[2] = state.waer[2];
        state.scase = "C2";
        cal_cc2(input, state);
        cal_cnh3p(state);
    }
}

void Solver::isrp2r(const Input& input, State& state) {
    state.clear_errors();
    state.scase = "2R";
    state.actmod = 2;

    bool tryliq = true;

    for (int it = 0; it < 2; ++it) {
        double sulratw = 2.0;
        if (tryliq && state.rh >= state.drnh4no3) {
            sulratw = getasr(state.waer[1], state.rh);
        }
        double sulrat = state.waer[2] / state.waer[1];

        if (sulratw <= sulrat) {
            state.scase = "N3";
            cal_n3(input, state);
        } else if (sulrat >= 1.0) {
            state.w[1] = state.waer[1];
            state.w[2] = state.waer[2];
            state.w[3] = state.waer[3];
            state.scase = "B4";
            cal_cb4(input, state);
            state.molal[3] = state.waer[3];
            state.molal[1] += state.waer[3];
            cal_nap(state);
            cal_cnh3p(state);
        } else {
            state.w[1] = state.waer[1];
            state.w[2] = state.waer[2];
            state.w[3] = state.waer[3];
            state.scase = "C2";
            cal_cc2(input, state);
            state.molal[3] = state.waer[3];
            state.molal[1] += state.waer[3];
            cal_nap(state);
            cal_cnh3p(state);
        }

        if (sulratw <= sulrat && sulrat < 2.0 && state.water <= state.tiny) {
            tryliq = false;
        } else {
            break;
        }
    }
}

//=======================================================================
// CASE S2 SPECIATION (NH4-SO4 deliquesced reverse analytical solver)
//=======================================================================
void Solver::cal_s2(const Input& input, State& state) {
    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    state.molalr[3] = std::min(state.waer[1], 0.5 * state.waer[2]); // (NH4)2SO4
    state.water     = state.molalr[3] / state.m0[3];
    state.water     = std::max(state.water, state.tiny);

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        state.a2  = state.xk2 * state.r * state.temp / (state.xkw * state.rh) * std::pow(state.gama[7] / state.gama[8], 2.0);
        double akw = state.xkw * state.rh * state.water * state.water;

        double nh4i  = state.waer[2];
        double so4i  = state.waer[1];
        double hso4i = 0.0;

        double hi = 0.0, ohi = 0.0;
        cal_cph(2.0 * so4i - nh4i, hi, ohi, state);

        double nh3aq = 0.0;
        if (hi < ohi) {
            double del = 0.0;
            cal_claq(nh4i, ohi, del, state); // Replaces CALCAMAQ with CALCLAQ Equivalent
            nh4i  = std::max(nh4i - del, 0.0);
            ohi   = std::max(ohi - del, state.tiny);
            nh3aq = del;
            hi    = akw / ohi;
        }

        double del = 0.0;
        cal_chs4(hi, so4i, 0.0, del, state);
        so4i  -= del;
        hi    -= del;
        hso4i  = del;

        double nh3gi = nh4i / hi / state.a2;

        state.molal[1] = hi;
        state.molal[2] = nh4i;
        state.molal[5] = so4i;
        state.molal[6] = hso4i;
        state.coh      = ohi;
        state.gasaq[0] = nh3aq;
        state.gnh3     = nh3gi;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act2();
        } else {
            break;
        }
    }
}
double Solver::funcs2(double x, const Input& input, State& state) {
    return 0.0;
}

//=======================================================================
// CASE N3 SPECIATION (NH4-SO4-NO3 deliquesced reverse analytical solver)
//=======================================================================
void Solver::cal_n3(const Input& input, State& state) {
    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    state.molalr[3] = std::min(state.waer[1], 0.5 * state.waer[2]); // (NH4)2SO4
    double aml5      = std::max(state.waer[2] - 2.0 * state.molalr[3], 0.0);
    state.molalr[4] = std::max(std::min(aml5, state.waer[3]), 0.0); // NH4NO3
    state.water     = state.molalr[3] / state.m0[3] + state.molalr[4] / state.m0[4];
    state.water     = std::max(state.water, state.tiny);

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        state.a2  = state.xk2 * state.r * state.temp / (state.xkw * state.rh) * std::pow(state.gama[7] / state.gama[8], 2.0);
        state.a3  = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
        double akw = state.xkw * state.rh * state.water * state.water;

        double nh4i  = state.waer[2];
        double no3i  = state.waer[3];
        double so4i  = state.waer[1];
        double hso4i = 0.0;

        double hi = 0.0, ohi = 0.0;
        cal_cph(2.0 * so4i + no3i - nh4i, hi, ohi, state);

        double nh3aq = 0.0;
        double no3aq = 0.0;
        double gg    = 2.0 * so4i + no3i - nh4i;

        if (hi < ohi) {
            double del = 0.0;
            cal_claq(-gg, nh4i, del, state); // Replaces CALCAMAQ2 with CALCLAQ Equivalent
            nh3aq = del;
            hi    = akw / ohi;
        } else {
            hi = 0.0;
            double del = 0.0;
            cal_niaq(gg, no3i, del, state); // Replaces CALCNIAQ2 with CALCNIAQ Equivalent
            no3aq = del;

            double del_hs = 0.0;
            cal_chs4(hi, so4i, 0.0, del_hs, state);
            so4i  -= del_hs;
            hi    -= del_hs;
            hso4i  = del_hs;
            ohi    = akw / hi;
        }

        state.molal[1] = hi;
        state.molal[2] = nh4i;
        state.molal[5] = so4i;
        state.molal[6] = hso4i;
        state.molal[3] = no3i;
        state.coh      = ohi;

        state.cnh42s4 = 0.0;
        state.cnh4no3 = 0.0;

        state.gasaq[0] = nh3aq;
        state.gasaq[2] = no3aq;

        state.ghno3 = hi * no3i / state.a3;
        state.gnh3  = nh4i / hi / state.a2;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act2();
        } else {
            break;
        }
    }
}
double Solver::funcn3(double x, const Input& input, State& state) {
    return 0.0;
}

void Solver::cal_nap(State& state) {
    if (state.water <= state.tiny) return;
    double delt = 0.0;
    cal_niaq(state.molal[3], state.molal[1], delt, state);
    state.a4 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
    state.gasaq[2] = delt;
    state.molal[1] -= delt;
    state.molal[3] -= delt;
    state.ghno3     = state.molal[1] * state.molal[3] / state.a4;
}

void Solver::isrp3r(const Input& input, State& state) {
    state.clear_errors();
    state.scase = "3R";
    state.actmod = 3;

    bool tryliq = true;

    for (int it = 0; it < 2; ++it) {
        double sulratw = 2.0;
        if (tryliq && state.rh >= state.drnh4no3) {
            double frso4 = state.waer[1] - state.waer[0] / 2.0;
            frso4 = std::max(frso4, state.tiny);
            double sri = getasr(frso4, state.rh);
            sulratw = (state.waer[0] + frso4 * sri) / state.waer[1];
            sulratw = std::min(sulratw, 2.0);
        }

        double sulrat = (state.waer[0] + state.waer[2]) / state.waer[1];
        double sodrat = state.waer[0] / state.waer[1];

        if (sulratw <= sulrat && sodrat < 2.0) {
            state.scase = "Q5";
            cal_q5(input, state);
        } else if (sulrat >= sulratw && sodrat >= 2.0) {
            state.scase = "R6";
            cal_r6(input, state);
        } else if (sulrat >= 1.0) {
            for (size_t i = 0; i < 8; ++i) state.w[i] = state.waer[i];
            state.scase = "I6";
            cal_ci6(input, state);
            cal_cnhp(input, state);
            cal_cnh3p(state);
        } else {
            for (size_t i = 0; i < 8; ++i) state.w[i] = state.waer[i];
            state.scase = "J3";
            cal_cj3(input, state);
            cal_cnhp(input, state);
            cal_cnh3p(state);
        }

        if (sulratw <= sulrat && sulrat < 2.0 && state.water <= state.tiny) {
            tryliq = false;
        } else {
            break;
        }
    }
}

//=======================================================================
// CASE Q5 SPECIATION (Marine deliquesced reverse analytical solver)
//=======================================================================
void Solver::cal_q5(const Input& input, State& state) {
    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    cal_q1a(input, state);

    state.psi1 = state.cna2so4;
    state.psi4 = state.cnh4cl;
    state.psi5 = state.cnh4no3;
    state.psi6 = state.cnh42s4;

    state.cal_cmr();

    double nh3aq = 0.0;
    double no3aq = 0.0;
    double claq  = 0.0;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double akw = state.xkw * state.rh * state.water * state.water;

        double nai  = state.waer[0];
        double so4i  = state.waer[1];
        double nh4i  = state.waer[2];
        double no3i  = state.waer[3];
        double cli   = state.waer[4];
        double hso4i = 0.0;

        double gg = 2.0 * so4i + no3i + cli - nai - nh4i;
        double hi = 0.0, ohi = 0.0;

        if (gg > state.tiny) {
            double bb = -gg;
            double cc = -akw;
            double dd = bb * bb - 4.0 * cc;
            hi  = 0.5 * (-bb + std::sqrt(dd));
            ohi = akw / hi;
        } else {
            double bb = gg;
            double cc = -akw;
            double dd = bb * bb - 4.0 * cc;
            ohi = 0.5 * (-bb + std::sqrt(dd));
            hi  = akw / ohi;
        }

        if (hi < ohi) {
            double del = 0.0;
            cal_claq(-gg, nh4i, del, state); // Replaces CALCAMAQ2 with CALCLAQ Equivalent
            nh3aq = del;
            hi    = akw / ohi;
            hso4i = 0.0;
        } else {
            double ggno3 = std::max(2.0 * so4i + no3i - nai - nh4i, 0.0);
            double ggcl  = std::max(gg - ggno3, 0.0);
            if (ggcl > state.tiny) {
                double del = 0.0;
                cal_claq(ggcl, cli, hi, state); // Replaces CALCCLAQ2 with CALCLAQ Equivalent (cal_claq takes double& delt)
                // Wait! Let's pass parameters correctly to CALCCLAQ2 (reproduced as Solver::cal_claq(cli, hi, delt, state))
                cal_claq(cli, hi, del, state);
                claq = del;
            }
            if (ggno3 > state.tiny) {
                if (ggcl <= state.tiny) hi = 0.0;
                double del = 0.0;
                cal_niaq(no3i, hi, del, state); // Replaces CALCNIAQ2 with CALCNIAQ equivalent
                no3aq = del;
            }

            double del = 0.0;
            cal_chs4(hi, so4i, 0.0, del, state);
            so4i  -= del;
            hi    -= del;
            hso4i  = del;
            ohi    = akw / hi;
        }

        state.molal[0] = nai;
        state.molal[1] = hi;
        state.molal[2] = nh4i;
        state.molal[3] = no3i;
        state.molal[4] = cli;
        state.molal[5] = so4i;
        state.molal[6] = hso4i;
        state.coh      = ohi;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act2();
        } else {
            break;
        }
    }

    state.a2 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
    state.a3 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
    state.a4 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);

    state.gnh3  = state.molal[2] / state.molal[1] / state.a2;
    state.ghno3 = state.molal[1] * state.molal[3] / state.a3;
    state.ghcl  = state.molal[1] * state.molal[4] / state.a4;

    state.gasaq[0] = nh3aq;
    state.gasaq[1] = claq;
    state.gasaq[2] = no3aq;

    state.cnh42s4 = 0.0;
    state.cnh4no3 = 0.0;
    state.cnh4cl  = 0.0;
    state.cnacl   = 0.0;
    state.cnano3  = 0.0;
    state.cna2so4 = 0.0;
}
double Solver::funcq5(double x, const Input& input, State& state) {
    return 0.0;
}

void Solver::cal_q1a(const Input& input, State& state) {
    state.cna2so4 = 0.5 * state.waer[0];
    double frso4   = std::max(state.waer[1] - state.cna2so4, 0.0);

    state.cnh42s4 = std::max(std::min(frso4, 0.5 * state.waer[2]), state.tiny);
    double frnh3   = std::max(state.waer[2] - 2.0 * state.cnh42s4, 0.0);

    state.cnh4no3 = std::min(frnh3, state.waer[3]);
    frnh3         = std::max(frnh3 - state.cnh4no3, 0.0);

    state.cnh4cl  = std::min(frnh3, state.waer[4]);
    frnh3         = std::max(frnh3 - state.cnh4cl, 0.0);

    state.water = 0.0;
    state.gnh3  = 0.0;
    state.ghno3 = 0.0;
    state.ghcl  = 0.0;
}

void Solver::isrp4r(const Input& input, State& state) {
    state.clear_errors();
    state.scase = "4R";
    state.actmod = 4;

    bool tryliq = true;

    for (int it = 0; it < 2; ++it) {
        double sulratw = 2.0;
        if (tryliq) {
            double frso4 = state.waer[1] - state.waer[0] / 2.0 - state.waer[5] - state.waer[6] / 2.0 - state.waer[7];
            frso4 = std::max(frso4, state.tiny);
            double sri = getasr(frso4, state.rh);
            sulratw = (state.waer[0] + frso4 * sri + state.waer[5] + state.waer[6] + state.waer[7]) / state.waer[1];
            sulratw = std::min(sulratw, 2.0);
        }

        double sulrat  = (state.waer[0] + state.waer[2] + state.waer[5] + state.waer[6] + state.waer[7]) / state.waer[1];
        double crnarat = (state.waer[0] + state.waer[5] + state.waer[6] + state.waer[7]) / state.waer[1];
        double crrat   = (state.waer[5] + state.waer[6] + state.waer[7]) / state.waer[1];

        if (sulrat >= sulratw) {
            if (crnarat < 2.0) {
                state.scase = "V7";
                cal_v7(input, state);
            } else if (crnarat >= 2.0 && crrat < 2.0) {
                state.scase = "U8";
                cal_u8(input, state);
            } else {
                state.scase = "W13";
                cal_w13(input, state);
            }
        } else if (sulrat >= 1.0) {
            for (size_t i = 0; i < 8; ++i) state.w[i] = state.waer[i];
            state.scase = "L9";
            cal_cl9(input, state);
            cal_cnhp(input, state);
            cal_cnh3p(state);
        } else {
            for (size_t i = 0; i < 8; ++i) state.w[i] = state.waer[i];
            state.scase = "K4";
            cal_ck4(input, state);
            cal_cnhp(input, state);
            cal_cnh3p(state);
        }

        if (sulrat >= sulratw && sulrat < 2.0 && state.water <= state.tiny) {
            tryliq = false;
        } else {
            break;
        }
    }
}

//=======================================================================
// CASE R6 SPECIATION (Sulfate poor, Sodium rich deliquesced reverse solver)
//=======================================================================
void Solver::cal_r6(const Input& input, State& state) {
    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    cal_r1a(input, state);

    state.psi1 = state.cna2so4;
    state.psi2 = state.cnano3;
    state.psi3 = state.cnacl;
    state.psi4 = state.cnh4cl;
    state.psi5 = state.cnh4no3;

    state.cal_cmr();

    double nh3aq = 0.0;
    double no3aq = 0.0;
    double claq  = 0.0;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double akw = state.xkw * state.rh * state.water * state.water;

        double nai  = state.waer[0];
        double so4i  = state.waer[1];
        double nh4i  = state.waer[2];
        double no3i  = state.waer[3];
        double cli   = state.waer[4];
        double hso4i = 0.0;

        double gg = 2.0 * so4i + no3i + cli - nai - nh4i;
        double hi = 0.0, ohi = 0.0;

        if (gg > state.tiny) {
            double bb = -gg;
            double cc = -akw;
            double dd = bb * bb - 4.0 * cc;
            hi  = 0.5 * (-bb + std::sqrt(dd));
            ohi = akw / hi;
        } else {
            double bb = gg;
            double cc = -akw;
            double dd = bb * bb - 4.0 * cc;
            ohi = 0.5 * (-bb + std::sqrt(dd));
            hi  = akw / ohi;
        }

        if (hi < ohi) {
            double del = 0.0;
            cal_claq(-gg, nh4i, del, state);
            nh3aq = del;
            hi    = akw / ohi;
        } else {
            double ggno3 = std::max(2.0 * so4i + no3i - nai - nh4i, 0.0);
            double ggcl  = std::max(gg - ggno3, 0.0);
            if (ggcl > state.tiny) {
                double del = 0.0;
                cal_claq(ggcl, cli, hi, state);
                claq = del;
            }
            if (ggno3 > state.tiny) {
                if (ggcl <= state.tiny) hi = 0.0;
                double del = 0.0;
                cal_niaq(ggno3, no3i, hi, state);
                no3aq = del;
            }

            double del = 0.0;
            cal_chs4(hi, so4i, 0.0, del, state);
            so4i  -= del;
            hi    -= del;
            hso4i  = del;
            ohi    = akw / hi;
        }

        state.molal[0] = nai;
        state.molal[1] = hi;
        state.molal[2] = nh4i;
        state.molal[3] = no3i;
        state.molal[4] = cli;
        state.molal[5] = so4i;
        state.molal[6] = hso4i;
        state.coh      = ohi;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act2();
        } else {
            break;
        }
    }

    state.a2 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
    state.a3 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
    state.a4 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);

    state.gnh3  = state.molal[2] / state.molal[1] / state.a2;
    state.ghno3 = state.molal[1] * state.molal[3] / state.a3;
    state.ghcl  = state.molal[1] * state.molal[4] / state.a4;

    state.gasaq[0] = nh3aq;
    state.gasaq[1] = claq;
    state.gasaq[2] = no3aq;

    state.cnh42s4 = 0.0;
    state.cnh4no3 = 0.0;
    state.cnh4cl  = 0.0;
    state.cnacl   = 0.0;
    state.cnano3  = 0.0;
    state.cna2so4 = 0.0;
}
double Solver::funcr6(double x, const Input& input, State& state) {
    return 0.0;
}

void Solver::cal_r1a(const Input& input, State& state) {
    state.cna2so4 = state.waer[1];
    double frna   = std::max(state.waer[0] - 2.0 * state.cna2so4, 0.0);

    state.cnh42s4 = 0.0;

    state.cnano3  = std::min(frna, state.waer[3]);
    double frno3  = std::max(state.waer[3] - state.cnano3, 0.0);
    frna          = std::max(frna - state.cnano3, 0.0);

    state.cnacl   = std::min(frna, state.waer[4]);
    double frcl   = std::max(state.waer[4] - state.cnacl, 0.0);
    frna          = std::max(frna - state.cnacl, 0.0);

    state.cnh4no3 = std::min(frno3, state.waer[2]);
    frno3         = std::max(frno3 - state.cnh4no3, 0.0);
    double frnh3  = std::max(state.waer[2] - state.cnh4no3, 0.0);

    state.cnh4cl  = std::min(frcl, frnh3);
    frcl          = std::max(frcl - state.cnh4cl, 0.0);
    frnh3         = std::max(frnh3 - state.cnh4cl, 0.0);

    state.water = 0.0;
    state.gnh3  = 0.0;
    state.ghno3 = 0.0;
    state.ghcl  = 0.0;
}

//=======================================================================
// CASE V7 SPECIATION (Sulfate poor, Cr+Na poor deliquesced reverse solver)
//=======================================================================
void Solver::cal_v7(const Input& input, State& state) {
    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    cal_v1a(input, state);

    state.psi1 = state.cna2so4;
    state.psi4 = state.cnh4cl;
    state.psi5 = state.cnh4no3;
    state.psi6 = state.cnh42s4;
    state.psi7 = state.ck2so4;
    state.psi8 = state.cmgso4;
    state.psi9 = state.ccaso4;

    state.cal_cmr();

    double nh3aq = 0.0;
    double no3aq = 0.0;
    double claq  = 0.0;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double akw = state.xkw * state.rh * state.water * state.water;

        double nai   = state.waer[0];
        double so4i  = std::max(state.waer[1] - state.waer[5], 0.0);
        double nh4i  = state.waer[2];
        double no3i  = state.waer[3];
        double cli   = state.waer[4];
        double cai   = 0.0;
        double ki    = state.waer[6];
        double mgi   = state.waer[7];
        double hso4i = 0.0;

        double gg = 2.0 * so4i + no3i + cli - nai - nh4i - 2.0 * cai - ki - 2.0 * mgi;
        double hi = 0.0, ohi = 0.0;

        if (gg > state.tiny) {
            double bb = -gg;
            double cc = -akw;
            double dd = bb * bb - 4.0 * cc;
            hi  = 0.5 * (-bb + std::sqrt(dd));
            ohi = akw / hi;
        } else {
            double bb = gg;
            double cc = -akw;
            double dd = bb * bb - 4.0 * cc;
            ohi = 0.5 * (-bb + std::sqrt(dd));
            hi  = akw / ohi;
        }

        double del = 0.0;
        if (hi > ohi) {
            cal_chs4(hi, so4i, 0.0, del, state);
        }

        so4i  -= del;
        hi    -= del;
        hso4i  = del;

        if (hi <= state.tiny) {
            hi  = std::sqrt(akw);
            ohi = akw / hi;
        } else {
            ohi = akw / hi;
        }

        state.molal[0] = nai;
        state.molal[1] = hi;
        state.molal[2] = nh4i;
        state.molal[3] = no3i;
        state.molal[4] = cli;
        state.molal[5] = so4i;
        state.molal[6] = hso4i;
        state.molal[7] = cai;
        state.molal[8] = ki;
        state.molal[9] = mgi;
        state.coh      = ohi;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act2();
        } else {
            break;
        }
    }

    state.a2 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
    state.a3 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
    state.a4 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);

    state.gnh3  = state.molal[2] / state.molal[1] / state.a2;
    state.ghno3 = state.molal[1] * state.molal[3] / state.a3;
    state.ghcl  = state.molal[1] * state.molal[4] / state.a4;

    state.gasaq[0] = nh3aq;
    state.gasaq[1] = claq;
    state.gasaq[2] = no3aq;

    state.cnh42s4 = 0.0;
    state.cnh4no3 = 0.0;
    state.cnh4cl  = 0.0;
    state.cna2so4 = 0.0;
    state.cmgso4  = 0.0;
    state.ck2so4  = 0.0;
    state.ccaso4  = std::min(state.waer[5], state.waer[1]);
}
double Solver::funcv7(double x, const Input& input, State& state) { return 0.0; }

void Solver::cal_v1a(const Input& input, State& state) {
    state.ccaso4  = std::min(state.waer[5], state.waer[1]);
    double frso4   = std::max(state.waer[1] - state.ccaso4, 0.0);
    double cafr    = std::max(state.waer[5] - state.ccaso4, 0.0);
    state.ck2so4  = std::min(0.5 * state.waer[6], frso4);
    double frk     = std::max(state.waer[6] - 2.0 * state.ck2so4, 0.0);
    frso4          = std::max(frso4 - state.ck2so4, 0.0);
    state.cna2so4 = std::min(0.5 * state.waer[0], frso4);
    double frna    = std::max(state.waer[0] - 2.0 * state.cna2so4, 0.0);
    frso4          = std::max(frso4 - state.cna2so4, 0.0);
    state.cmgso4  = std::min(state.waer[7], frso4);
    double frmg    = std::max(state.waer[7] - state.cmgso4, 0.0);
    frso4          = std::max(frso4 - state.cmgso4, 0.0);
    state.cnh42s4 = std::max(std::min(frso4, 0.5 * state.waer[2]), state.tiny);
    double frnh3   = std::max(state.waer[2] - 2.0 * state.cnh42s4, 0.0);

    state.cnh4no3 = std::min(frnh3, state.waer[3]);
    frnh3         = std::max(frnh3 - state.cnh4no3, 0.0);

    state.cnh4cl  = std::min(frnh3, state.waer[4]);
    frnh3         = std::max(frnh3 - state.cnh4cl, 0.0);

    state.water = 0.0;
    state.gnh3  = 0.0;
    state.ghno3 = 0.0;
    state.ghcl  = 0.0;
}

//=======================================================================
// CASE U8 SPECIATION (Sulfate poor, Crustal+Sodium rich reverse solver)
//=======================================================================
void Solver::cal_u8(const Input& input, State& state) {
    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    cal_u1a(input, state);

    state.psi1 = state.cna2so4;
    state.psi2 = state.cnano3;
    state.psi3 = state.cnacl;
    state.psi4 = state.cnh4cl;
    state.psi5 = state.cnh4no3;
    state.psi7 = state.ck2so4;
    state.psi8 = state.cmgso4;
    state.psi9 = state.ccaso4;

    state.cal_cmr();

    double nh3aq = 0.0;
    double no3aq = 0.0;
    double claq  = 0.0;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double akw = state.xkw * state.rh * state.water * state.water;

        double nai   = state.waer[0];
        double so4i  = std::max(state.waer[1] - state.waer[5], 0.0);
        double nh4i  = state.waer[2];
        double no3i  = state.waer[3];
        double cli   = state.waer[4];
        double cai   = 0.0;
        double ki    = state.waer[6];
        double mgi   = state.waer[7];
        double hso4i = 0.0;

        double gg = 2.0 * so4i + no3i + cli - nai - nh4i - 2.0 * cai - ki - 2.0 * mgi;
        double hi = 0.0, ohi = 0.0;

        if (gg > state.tiny) {
            double bb = -gg;
            double cc = -akw;
            double dd = bb * bb - 4.0 * cc;
            hi  = 0.5 * (-bb + std::sqrt(dd));
            ohi = akw / hi;
        } else {
            double bb = gg;
            double cc = -akw;
            double dd = bb * bb - 4.0 * cc;
            ohi = 0.5 * (-bb + std::sqrt(dd));
            hi  = akw / ohi;
        }

        if (hi <= state.tiny) hi = std::sqrt(akw);

        double del = 0.0;
        if (hi > ohi) {
            cal_chs4(hi, so4i, 0.0, del, state);
        }

        so4i  -= del;
        hi    -= del;
        hso4i  = del;
        ohi    = akw / hi;

        if (hi <= state.tiny) {
            hi  = std::sqrt(akw);
            ohi = akw / hi;
        }

        state.molal[0] = nai;
        state.molal[1] = hi;
        state.molal[2] = nh4i;
        state.molal[3] = no3i;
        state.molal[4] = cli;
        state.molal[5] = so4i;
        state.molal[6] = hso4i;
        state.molal[7] = cai;
        state.molal[8] = ki;
        state.molal[9] = mgi;
        state.coh      = ohi;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act2();
        } else {
            break;
        }
    }

    state.a2 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
    state.a3 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
    state.a4 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);

    state.gnh3  = state.molal[2] / state.molal[1] / state.a2;
    state.ghno3 = state.molal[1] * state.molal[3] / state.a3;
    state.ghcl  = state.molal[1] * state.molal[4] / state.a4;

    state.gasaq[0] = nh3aq;
    state.gasaq[1] = claq;
    state.gasaq[2] = no3aq;

    state.cnh42s4 = 0.0;
    state.cnh4no3 = 0.0;
    state.cnh4cl  = 0.0;
    state.cnacl   = 0.0;
    state.cnano3  = 0.0;
    state.cna2so4 = 0.0;
    state.cmgso4  = 0.0;
    state.ck2so4  = 0.0;
    state.ccaso4  = std::min(state.waer[5], state.waer[1]);
}
double Solver::funcu8(double x, const Input& input, State& state) { return 0.0; }

void Solver::cal_u1a(const Input& input, State& state) {
    state.ccaso4  = std::min(state.waer[5], state.waer[1]);
    double frso4   = std::max(state.waer[1] - state.ccaso4, 0.0);
    double cafr    = std::max(state.waer[5] - state.ccaso4, 0.0);
    state.ck2so4  = std::min(0.5 * state.waer[6], frso4);
    double frk     = std::max(state.waer[6] - 2.0 * state.ck2so4, 0.0);
    frso4          = std::max(frso4 - state.ck2so4, 0.0);
    state.cmgso4  = std::min(state.waer[7], frso4);
    double frmg    = std::max(state.waer[7] - state.cmgso4, 0.0);
    frso4          = std::max(frso4 - state.cmgso4, 0.0);
    state.cna2so4 = std::max(frso4, 0.0);
    double frna    = std::max(state.waer[0] - 2.0 * state.cna2so4, 0.0);

    state.cnh42s4 = 0.0;

    state.cnano3  = std::min(frna, state.waer[3]);
    double frno3  = std::max(state.waer[3] - state.cnano3, 0.0);
    frna          = std::max(frna - state.cnano3, 0.0);

    state.cnacl   = std::min(frna, state.waer[4]);
    double frcl   = std::max(state.waer[4] - state.cnacl, 0.0);
    frna          = std::max(frna - state.cnacl, 0.0);

    state.cnh4no3 = std::min(frno3, state.waer[2]);
    frno3         = std::max(frno3 - state.cnh4no3, 0.0);
    double frnh3  = std::max(state.waer[2] - state.cnh4no3, 0.0);

    state.cnh4cl  = std::min(frcl, frnh3);
    frcl          = std::max(frcl - state.cnh4cl, 0.0);
    frnh3         = std::max(frnh3 - state.cnh4cl, 0.0);

    state.water = 0.0;
    state.gnh3  = 0.0;
    state.ghno3 = 0.0;
    state.ghcl  = 0.0;
}

//=======================================================================
// CASE W13 SPECIATION (Sulfate poor, Crustal rich deliquesced reverse solver)
//=======================================================================
void Solver::cal_w13(const Input& input, State& state) {
    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    cal_w1a(input, state);

    state.psi1  = state.cna2so4;
    state.psi5  = state.cnh4cl;
    state.psi6  = state.cnh4no3;
    state.psi7  = state.cnacl;
    state.psi8  = state.cnano3;
    state.psi9  = state.ck2so4;
    state.psi10 = state.cmgso4;
    state.psi11 = state.ccaso4;
    state.psi12 = state.ccano32;
    state.psi13 = state.ckno3;
    state.psi14 = state.ckcl;
    state.psi15 = state.cmgno32;
    state.psi16 = state.cmgcl2;
    state.psi17 = state.ccacl2;

    state.cal_cmr();

    double nh3aq = 0.0;
    double no3aq = 0.0;
    double claq  = 0.0;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double akw = state.xkw * state.rh * state.water * state.water;

        double nai   = state.waer[0];
        double so4i  = std::max(state.waer[1] - state.waer[5], 0.0);
        double nh4i  = state.waer[2];
        double no3i  = state.waer[3];
        double cli   = state.waer[4];
        double cai   = 0.0;
        double ki    = state.waer[6];
        double mgi   = state.waer[7];
        double hso4i = 0.0;

        double gg = 2.0 * so4i + no3i + cli - nai - nh4i - 2.0 * cai - ki - 2.0 * mgi;
        double hi = 0.0, ohi = 0.0;

        if (gg > state.tiny) {
            double bb = -gg;
            double cc = -akw;
            double dd = bb * bb - 4.0 * cc;
            hi  = 0.5 * (-bb + std::sqrt(dd));
            ohi = akw / hi;
        } else {
            double bb = gg;
            double cc = -akw;
            double dd = bb * bb - 4.0 * cc;
            ohi = 0.5 * (-bb + std::sqrt(dd));
            hi  = akw / ohi;
        }

        double del = 0.0;
        if (hi > ohi) {
            cal_chs4(hi, so4i, 0.0, del, state);
        }

        so4i  -= del;
        hi    -= del;
        hso4i  = del;
        ohi    = akw / hi;

        if (hi <= state.tiny) {
            hi  = std::sqrt(akw);
            ohi = akw / hi;
        }

        state.molal[0] = nai;
        state.molal[1] = hi;
        state.molal[2] = nh4i;
        state.molal[3] = no3i;
        state.molal[4] = cli;
        state.molal[5] = so4i;
        state.molal[6] = hso4i;
        state.molal[7] = cai;
        state.molal[8] = ki;
        state.molal[9] = mgi;
        state.coh      = ohi;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act2();
        } else {
            break;
        }
    }

    state.a2 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
    state.a3 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
    state.a4 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);

    state.gnh3  = state.molal[2] / state.molal[1] / state.a2;
    state.ghno3 = state.molal[1] * state.molal[3] / state.a3;
    state.ghcl  = state.molal[1] * state.molal[4] / state.a4;

    state.gasaq[0] = nh3aq;
    state.gasaq[1] = claq;
    state.gasaq[2] = no3aq;

    state.cnh42s4 = 0.0;
    state.cnh4no3 = 0.0;
    state.cnh4cl  = 0.0;
    state.cnacl   = 0.0;
    state.cnano3  = 0.0;
    state.cmgso4  = 0.0;
    state.ck2so4  = 0.0;
    state.ccaso4  = std::min(state.waer[5], state.waer[1]);
    state.ccano32 = 0.0;
    state.ckno3   = 0.0;
    state.ckcl    = 0.0;
    state.cmgno32 = 0.0;
    state.cmgcl2  = 0.0;
    state.ccacl2  = 0.0;
}
double Solver::funcw13(double x, const Input& input, State& state) { return 0.0; }

void Solver::cal_w1a(const Input& input, State& state) {
    state.ccaso4  = std::min(state.waer[1], state.waer[5]);
    double cafr    = std::max(state.waer[5] - state.ccaso4, 0.0);
    double so4fr   = std::max(state.waer[1] - state.ccaso4, 0.0);
    state.ck2so4  = std::min(so4fr, 0.5 * state.waer[6]);
    double frk     = std::max(state.waer[6] - 2.0 * state.ck2so4, 0.0);
    so4fr          = std::max(so4fr - state.ck2so4, 0.0);
    state.cmgso4  = so4fr;
    double frmg    = std::max(state.waer[7] - state.cmgso4, 0.0);
    state.cnacl   = std::min(state.waer[0], state.waer[4]);
    double frna    = std::max(state.waer[0] - state.cnacl, 0.0);
    double clfr    = std::max(state.waer[4] - state.cnacl, 0.0);
    state.ccacl2  = std::min(cafr, 0.5 * clfr);
    cafr           = std::max(cafr - state.ccacl2, 0.0);
    clfr           = std::max(state.waer[4] - 2.0 * state.ccacl2, 0.0);
    state.ccano32 = std::min(cafr, 0.5 * state.waer[3]);
    cafr           = std::max(cafr - state.ccano32, 0.0);
    double frno3   = std::max(state.waer[3] - 2.0 * state.ccano32, 0.0);
    state.cmgcl2  = std::min(frmg, 0.5 * clfr);
    frmg           = std::max(frmg - state.cmgcl2, 0.0);
    clfr           = std::max(clfr - 2.0 * state.cmgcl2, 0.0);
    state.cmgno32 = std::min(frmg, 0.5 * frno3);
    frmg           = std::max(frmg - state.cmgno32, 0.0);
    frno3          = std::max(frno3 - 2.0 * state.cmgno32, 0.0);
    state.cnano3  = std::min(frna, frno3);
    frna           = std::max(frna - state.cnano3, 0.0);
    frno3          = std::max(frno3 - state.cnano3, 0.0);
    state.ckcl    = std::min(frk, clfr);
    frk            = std::max(frk - state.ckcl, 0.0);
    clfr           = std::max(clfr - state.ckcl, 0.0);
    state.ckno3   = std::min(frk, frno3);
    frk            = std::max(frk - state.ckno3, 0.0);
    frno3          = std::max(frno3 - state.ckno3, 0.0);

    state.water = 0.0;
    state.gnh3  = 0.0;
    state.ghno3 = 0.0;
    state.ghcl  = 0.0;
}

} // namespace Isorropia
