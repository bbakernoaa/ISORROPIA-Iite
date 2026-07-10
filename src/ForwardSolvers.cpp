#include "Isorropia/Solver.hpp"
#include <cmath>
#include <algorithm>

namespace Isorropia {

void Solver::isrp1f(const Input& input, State& state) {
    state.clear_errors();
    state.actmod = 1; // Pre-set standard active activity coefficients model
    
    // SULRAT = W(3) / W(2) (Total Ammonia / Total Sulfate)
    double sulrat = state.w[2] / state.w[1];

    if (sulrat >= 2.0) {
        // Sulfate Poor Range
        double dc = state.w[2] - 2.001 * state.w[1];
        state.w[2] += std::max(-dc, 0.0);
        
        state.scase = "A2";
        cal_ca2(input, state);
    } 
    else if (sulrat >= 1.0) {
        // Sulfate Rich (No Free Acid) Range
        state.scase = "B4";
        cal_cb4(input, state);
        cal_cnh3(input, state);
    } 
    else {
        // Sulfate Rich (Free Acid) Range
        state.scase = "C2";
        cal_cc2(input, state);
        cal_cnh3(input, state);
    }
}

//=======================================================================
// CASE A2 - SULFATE POOR METASTABLE LIQUID Speciation
//=======================================================================
void Solver::cal_ca2(const Input& input, State& state) {
    state.calaou = true;
    double omelo = state.tiny;
    double omehi = 2.0 * state.w[1]; // 2.0 * Total Sulfate

    // 1. Initialize water content
    state.molal[5] = state.w[1]; // SO4-- is initially set to total sulfate
    state.molal[6] = 0.0;        // HSO4- is initially zero
    state.cal_cmr();             // Recalculate water content

    // 2. Initial values for Bisection
    double x1 = omehi;
    double eps = 1e-6;
    double y1 = funca2(x1, input, state);
    if (std::abs(y1) <= eps) return;

    // 3. Root Tracking across intervals
    double x2 = x1;
    double y2 = y1;
    int ndiv = 5;
    double dx = (omehi - omelo) / static_cast<double>(ndiv);

    bool sign_changed = false;
    for (int i = 1; i <= ndiv; ++i) {
        x2 = std::max(x1 - dx, omelo);
        y2 = funca2(x2, input, state);
        if ((y1 < 0.0 && y2 > 0.0) || (y1 > 0.0 && y2 < 0.0)) {
            sign_changed = true;
            break;
        }
        x1 = x2;
        y1 = y2;
    }

    if (!sign_changed) {
        if (std::abs(y2) <= eps) {
            return;
        } else {
            state.push_error(1, "CALCA2: NO SOLUTION");
            return;
        }
    }

    // 4. Perform Bisection
    int maxit = 100;
    for (int i = 1; i <= maxit; ++i) {
        double x3 = 0.5 * (x1 + x2);
        state.rstgamp();
        double y3 = funca2(x3, input, state);

        if (((y1 < 0.0 && y3 <= 0.0) || (y1 > 0.0 && y3 >= 0.0)) == false) {
            // Sign change between y1 and y3
            y2 = y3;
            x2 = x3;
        } else {
            // Sign change between y3 and y2
            y1 = y3;
            x1 = x3;
        }

        if (std::abs(x2 - x1) <= eps * x1) {
            // Converged
            double final_x3 = 0.5 * (x1 + x2);
            state.rstgamp();
            funca2(final_x3, input, state);
            return;
        }
    }

    state.push_error(2, "CALCA2: NO CONVERGENCE");
    double final_x3 = 0.5 * (x1 + x2);
    state.rstgamp();
    funca2(final_x3, input, state);
}

double Solver::funca2(double omegi, const Input& input, State& state) {
    state.frst = true;
    state.calain = true;
    double psi = state.w[1]; // Initial amount of sulfate in solution

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; sweep++) {
        // A1 = XK1 * WATER / GAMA(7) * (GAMA(8) / GAMA(7))^2
        double a1 = state.xk1 * state.water / state.gama[6] * std::pow(state.gama[7] / state.gama[6], 2.0);
        // A2 = XK2 * R * TEMP / XKW * (GAMA(8) / GAMA(9))^2
        double a2 = state.xk2 * state.r * state.temp / state.xkw * std::pow(state.gama[7] / state.gama[8], 2.0);
        // A3 = XKW * RH * WATER * WATER
        double a3 = state.xkw * state.rh * state.water * state.water;

        double lamda = psi / (a1 / omegi + 1.0);
        double zeta = a3 / omegi;

        // Populate liquid species concentrations
        state.molal[1] = omegi; // H+ -> maps to MOLAL(2)
        state.molal[5] = std::max(psi - lamda, state.tiny); // SO4-- -> maps to MOLAL(6)
        state.molal[2] = std::max(state.w[2] / (1.0 / a2 / omegi + 1.0), 2.0 * state.molal[5]); // NH4+ -> maps to MOLAL(3)
        state.molal[6] = lamda; // HSO4- -> maps to MOLAL(7)

        state.gnh3 = std::max(state.w[2] - state.molal[2], state.tiny);
        state.coh  = zeta;

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act1();
        } else {
            break;
        }
    }

    double denom = (2.0 * state.molal[5] + state.molal[6]);
    return (state.molal[2] / denom - 1.0) + state.molal[1] / denom;
}

//=======================================================================
// CASE B4 - SULFATE RICH Speciation
//=======================================================================
void Solver::cal_cb4(const Input& input, State& state) {
    state.frst = true;
    state.calain = true;
    state.calaou = true;

    // 1. Dry material balance to initialize ZSR water estimation
    cal_cb1a(input, state);
    state.molalr[12] = state.clc;
    state.molalr[8]  = state.cnh4hs4;
    state.molalr[3]  = state.cnh42s4;
    state.clc     = 0.0;
    state.cnh4hs4 = 0.0;
    state.cnh42s4 = 0.0;

    state.water = state.molalr[12] / state.m0[12] + state.molalr[8] / state.m0[8] + state.molalr[3] / state.m0[3];
    state.molal[2] = state.w[2]; // NH4+ is total Ammonia

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double ak1 = state.xk1 * std::pow(state.gama[7] / state.gama[6], 2.0) * (state.water / state.gama[6]);
        double bet = state.w[1];
        double gam = state.molal[2];

        double bb = bet + ak1 - gam;
        double cc = -ak1 * bet;
        double dd = bb * bb - 4.0 * cc;

        // Speciations calculations
        state.molal[5] = std::max(state.tiny, std::min(0.5 * (-bb + std::sqrt(dd)), state.w[1])); // SO4--
        state.molal[6] = std::max(state.tiny, std::min(state.w[1] - state.molal[5], state.w[1])); // HSO4-
        state.molal[1] = std::max(state.tiny, std::min(ak1 * state.molal[6] / state.molal[5], state.w[1])); // H+

        state.cal_cmr();

        if (!state.calain) break;
        if (state.actmod == 2) {
            state.cal_act2();
        } else {
            state.cal_act1();
        }
    }
}

void Solver::cal_cb1a(const Input& input, State& state) {
    double x = 2.0 * state.w[1] - state.w[2]; // Equivalent NH4HSO4
    double y = state.w[2] - state.w[1];       // Equivalent (NH4)2SO4

    if (x <= y) {
        state.clc     = x;
        state.cnh4hs4 = 0.0;
        state.cnh42s4 = y - x;
    } else {
        state.clc     = y;
        state.cnh4hs4 = x - y;
        state.cnh42s4 = 0.0;
    }
}

//=======================================================================
// CASE C2 - SULFATE RICH WITH FREE ACID Speciation
//=======================================================================
void Solver::cal_cc2(const Input& input, State& state) {
    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    double lamda = state.w[2]; // NH4HSO4 initially in solution
    double psi   = state.w[1] - state.w[2]; // Free H2SO4 in solution

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double parm = state.water * state.xk1 / state.gama[6] * std::pow(state.gama[7] / state.gama[6], 2.0);
        double bb   = psi + parm;
        double cc   = -parm * (lamda + psi);
        double kapa = 0.5 * (-bb + std::sqrt(bb * bb - 4.0 * cc));

        state.molal[1] = psi + kapa; // H+
        state.molal[2] = lamda;      // NH4+
        state.molal[5] = kapa;       // SO4--
        state.molal[6] = std::max(lamda + psi - kapa, state.tiny); // HSO4-

        state.ch2so4 = std::max(state.molal[5] + state.molal[6] - state.molal[2], 0.0); // Free H2SO4
        state.cal_cmr();

        if (!state.calain) break;
        if (state.actmod == 2) {
            state.cal_act2();
        } else {
            state.cal_act1();
        }
    }
}

//=======================================================================
// GAS-LIQUID AMMONIA BINDINGS (CALCNH3 Equivalent)
//=======================================================================
void Solver::cal_cnh3(const Input& input, State& state) {
    if (state.water <= state.tiny) return;

    // A1 = (XK2/XKW) * R * TEMP * (GAMA(10)/GAMA(5))^2
    // GAMA(10) -> state.gama[9] (HNO3aq), GAMA(5) -> state.gama[4] (NH4Cl)
    double a1 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
    double chi1 = state.molal[2]; // NH4+ (index 2)
    double chi2 = state.molal[1]; // H+ (index 1)

    double bb = chi2 + 1.0 / a1;
    double cc = -chi1 / a1;
    double diak = std::sqrt(bb * bb - 4.0 * cc);
    double psi = 0.5 * (-bb + diak);
    psi = std::max(state.tiny, std::min(psi, chi1));

    state.gnh3 = psi;
    state.molal[2] = chi1 - psi; // NH4+
    state.molal[1] = chi2 + psi; // H+
}

void Solver::isrp2f(const Input& input, State& state) {
    state.clear_errors();
    state.actmod = 2; // Pre-set standard active activity coefficients model
    double sulrat = state.w[2] / state.w[1];

    if (sulrat >= 2.0) {
        state.scase = "D3";
        cal_cd3(input, state);
    } 
    else if (sulrat >= 1.0) {
        state.scase = "E4"; // Set F77 active case prefix 'E' before speciation
        cal_cb4(input, state);
        cal_cna(input, state);
    } 
    else {
        state.scase = "F2"; // Set F77 active case prefix 'F' before speciation
        cal_cc2(input, state);
        cal_cna(input, state);
    }
}

//=======================================================================
// CASE D3 - SULFATE POOR METASTABLE Speciation
//=======================================================================
void Solver::cal_cd3(const Input& input, State& state) {
    // 1. Dry composition material balance
    cal_cd1a(input, state);

    // Save dry compositions for bisection
    state.chi1 = state.cnh4no3;
    state.chi2 = state.cnh42s4;
    state.chi3 = state.ghno3;
    state.chi4 = state.gnh3;

    state.psi1 = state.cnh4no3;
    state.psi2 = state.cnh42s4;
    state.psi3 = 0.0;
    state.psi4 = 0.0;

    // 2. Initial water estimation content
    state.molal[5] = state.psi2;               // SO4-- (index 5)
    state.molal[6] = 0.0;                      // HSO4- (index 6)
    state.molal[2] = state.psi1;               // NH4+ (index 2)
    state.molal[3] = state.psi1;               // NO3- (index 3)
    state.cal_cmr();

    state.calaou = true;
    double psi4lo = state.tiny;
    double psi4hi = state.chi4;
    double eps = 1e-6;
    int ndiv = 5;

    // Initial values for bisection
GOTO_60:
    double x1 = psi4lo;
    state.rstgamp();
    double y1 = funcd3(x1, input, state);
    if (std::abs(y1) <= eps) return;
    double ylo = y1;

    // Root Tracking across divisions
    double dx = (psi4hi - psi4lo) / static_cast<double>(ndiv);
    double x2 = x1;
    double y2 = y1;
    bool sign_changed = false;

    for (int i = 1; i <= ndiv; ++i) {
        x2 = x1 + dx;
        state.rstgamp();
        y2 = funcd3(x2, input, state);
        if (y1 < 0.0 && y2 > 0.0) {
            sign_changed = true;
            break;
        }
        x1 = x2;
        y1 = y2;
    }

    if (!sign_changed) {
        double yhi = y1;
        if (std::abs(y2) < eps) {
            return;
        }
        else if (ylo < 0.0 && yhi < 0.0) {
            double p4 = state.tiny;
            state.rstgamp();
            funcd3(p4, input, state);
            goto GOTO_50;
        }
        else if (ylo > 0.0 && yhi > 0.0) {
            psi4hi = psi4lo;
            psi4lo = psi4lo - 0.1 * (state.psi1 + state.psi2);
            if (psi4lo < -(state.psi1 + state.psi2)) {
                state.push_error(1, "CALCD3: NO SOLUTION");
                return;
            } else {
                state.molal[5] = state.psi2;
                state.molal[6] = 0.0;
                state.molal[2] = state.psi1;
                state.molal[3] = state.psi1;
                state.cal_cmr();
                goto GOTO_60; // Redo root tracking
            }
        }
    }

    // Perform Bisection
    {
        int maxit = 100;
        for (int i = 1; i <= maxit; ++i) {
            double x3 = 0.5 * (x1 + x2);
            state.rstgamp();
            double y3 = funcd3(x3, input, state);

            if (((y1 < 0.0 && y3 <= 0.0) || (y1 > 0.0 && y3 >= 0.0)) == false) {
                y2 = y3;
                x2 = x3;
            } else {
                y1 = y3;
                x1 = x3;
            }

            if (std::abs(x2 - x1) <= eps * std::abs(x1)) {
                double final_x3 = 0.5 * (x1 + x2);
                state.rstgamp();
                funcd3(final_x3, input, state);
                goto GOTO_50;
            }
        }
    }

    state.push_error(2, "CALCD3: NO CONVERGENCE");
    {
        double final_x3 = 0.5 * (x1 + x2);
        state.rstgamp();
        funcd3(final_x3, input, state);
    }

GOTO_50:
    if (state.molal[1] > state.tiny) {
        double delta = 0.0;
        cal_chs4(state.molal[1], state.molal[5], 0.0, delta, state);
        state.molal[1] -= delta;
        state.molal[5] -= delta;
        state.molal[6] = delta;
    }
}

double Solver::funcd3(double p4, const Input& input, State& state) {
    state.frst = true;
    state.calain = true;
    state.psi4 = p4;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        // A2 = XK7 * (WATER / GAMA(4))^3
        double a2 = state.xk7 * std::pow(state.water / state.gama[3], 3.0);
        
        // A3 = XK4 * R * TEMP * (WATER / GAMA(10))^2
        double a3 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
        
        // A4 = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
        double a4 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
        
        // A7 = XKW * RH * WATER * WATER
        double a7 = state.xkw * state.rh * state.water * state.water;

        // PSI3 calculations
        double psi3 = a3 * a4 * state.chi3 * (state.chi4 - state.psi4) - state.psi1 * (2.0 * state.psi2 + state.psi1 + state.psi4);
        psi3 = psi3 / (a3 * a4 * (state.chi4 - state.psi4) + 2.0 * state.psi2 + state.psi1 + state.psi4);
        psi3 = std::min(std::max(psi3, 0.0), state.chi3);
        state.psi3 = psi3;

        double bb = state.psi4 - state.psi3;
        double denm = bb + std::sqrt(bb * bb + 4.0 * a7);
        if (denm <= state.tiny) {
            double abb = std::abs(bb);
            denm = (bb + abb) + 2.0 * a7 / abb;
        }
        double ahi = 2.0 * a7 / denm;

        // Speciation population
        state.molal[1] = ahi;                                   // H+ (index 1)
        state.molal[2] = state.psi1 + state.psi4 + 2.0 * state.psi2; // NH4+ (index 2)
        state.molal[5] = state.psi2;                            // SO4-- (index 5)
        state.molal[6] = 0.0;                                   // HSO4- (index 6)
        state.molal[3] = state.psi3 + state.psi1;               // NO3- (index 3)

        state.cnh42s4 = state.chi2 - state.psi2;
        state.cnh4no3 = 0.0;
        state.ghno3   = state.chi3 - state.psi3;
        state.gnh3    = state.chi4 - state.psi4;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act2();
        } else {
            break;
        }
    }

    double a4_val = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
    return state.molal[2] / state.molal[1] / std::max(state.gnh3, state.tiny) / a4_val - 1.0;
}

void Solver::cal_cd1a(const Input& input, State& state) {
    double parm = state.xk10 / (state.r * state.temp) / (state.r * state.temp);
    state.cnh42s4 = state.w[1]; // Total Sulfate
    double x = std::max(0.0, std::min(state.w[2] - 2.0 * state.cnh42s4, state.w[3])); // Total HNO3
    double ps = std::max(state.w[2] - x - 2.0 * state.cnh42s4, 0.0);
    double om = std::max(state.w[3] - x, 0.0);

    double omps = om + ps;
    double diak = std::sqrt(omps * omps + 4.0 * parm);
    double ze = std::min(x, 0.5 * (-omps + diak));

    state.cnh4no3 = x - ze;
    state.gnh3 = ps + ze;
    state.ghno3 = om + ze;
}

void Solver::cal_chs4(double hi, double so4i, double hso4i, double& delta, State& state) {
    if (state.water <= 10.0 * state.tiny) {
        delta = 0.0;
        return;
    }

    double a8 = state.xk1 * state.water / state.gama[6] * std::pow(state.gama[7] / state.gama[6], 2.0);
    double bb = -(hi + so4i + a8);
    double cc = hi * so4i - hso4i * a8;
    double dd = bb * bb - 4.0 * cc;

    if (dd >= 0.0) {
        double sqdd = std::sqrt(dd);
        double delta1 = 0.5 * (-bb + sqdd);
        double delta2 = 0.5 * (-bb - sqdd);
        if (hso4i <= state.tiny) {
            delta = delta2;
        } else if (hi * so4i >= a8 * hso4i) {
            delta = delta2;
        } else {
            delta = delta1;
        }
    } else {
        delta = 0.0;
    }
}

void Solver::cal_cna(const Input& input, State& state) {
    double x = state.w[3]; // Total HNO3
    double delt = 0.0;
    if (state.water > state.tiny) {
        double kapa = state.molal[1]; // H+
        double alfa = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
        double diak = std::sqrt((kapa + alfa) * (kapa + alfa) + 4.0 * alfa * x);
        delt = 0.5 * (-(kapa + alfa) + diak);
    }

    state.ghno3 = std::max(x - delt, 0.0);
    state.molal[3] = delt; // NO3-
    state.molal[1] = state.molal[1] + delt; // H+
}

void Solver::cal_cha(const Input& input, State& state) {
    double x = state.w[4]; // Total HCl -> Component Cl index 4
    double delt = 0.0;
    if (state.water > state.tiny) {
        double kapa = state.molal[1]; // H+
        double alfa = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0); // GAMA(11) -> index 10
        double diak = std::sqrt((kapa + alfa) * (kapa + alfa) + 4.0 * alfa * x);
        delt = 0.5 * (-(kapa + alfa) + diak);
    }

    state.ghcl = std::max(x - delt, 0.0);
    state.molal[4] = delt; // Cl- -> index 4
    state.molal[1] = state.molal[1] + delt; // H+
}

void Solver::cal_cnha(const Input& input, State& state) {
    if (state.water <= state.tiny) {
        state.ghcl  = std::max(state.w[4] - state.molal[4], state.tiny);
        state.ghno3 = std::max(state.w[3] - state.molal[3], state.tiny);
        return;
    }

    if (state.w[4] <= state.tiny && state.w[3] <= state.tiny) {
        return;
    }
    else if (state.w[4] <= state.tiny) {
        cal_cna(input, state);
        return;
    }
    else if (state.w[3] <= state.tiny) {
        cal_cha(input, state);
        return;
    }

    double a3 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0); // GAMA(10)
    double a4 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0); // GAMA(11)

    double delcl = 0.0;
    double delno = 0.0;

    double omega = state.molal[1]; // H+
    double chi3  = state.w[3];     // HNO3
    double chi4  = state.w[4];     // HCl

    double c1 = a3 * chi3;
    double c2 = a4 * chi4;
    double c3 = a3 - a4;

    double m1 = (c1 + c2 + (omega + a4) * c3) / c3;
    double m2 = ((omega + a4) * c2 - a4 * c3 * chi4) / c3;
    double m3 = -a4 * c2 * chi4 / c3;

    int islv = 1;
    poly3(m1, m2, m3, delcl, islv, state);
    if (islv != 0) {
        delcl = state.tiny;
        state.push_error(22, "CALCNHA: NO SOLUTION FOR HCL DISS");
    }
    delcl = std::min(delcl, chi4);

    delno = c1 * delcl / (c2 + c3 * delcl);
    delno = std::min(delno, chi3);

    if (delcl < 0.0 || delno < 0.0 || delcl > chi4 || delno > chi3) {
        delcl = state.tiny;
        delno = state.tiny;
        state.push_error(22, "CALCNHA: SOLUTION BOUNDS WARN");
    }

    state.molal[1] += (delno + delcl); // H+
    state.molal[4] += delcl;           // Cl-
    state.molal[3] += delno;           // NO3-

    state.ghcl  = std::max(state.w[4] - state.molal[4], state.tiny);
    state.ghno3 = std::max(state.w[3] - state.molal[3], state.tiny);
}

void Solver::poly3(double a1, double a2, double a3, double& root, int& islv, State& state) {
    double eps = 1e-50;
    double expon = 1.0 / 3.0;
    double pi = 3.14159265358979323846;
    double thet1 = 120.0 / 180.0;
    double thet2 = 240.0 / 180.0;

    std::array<double, 3> x = {0.0};
    int ix = 1;

    if (std::abs(a3) <= eps) {
        islv = 1;
        ix = 1;
        x[0] = 0.0;
        double d = a1 * a1 - 4.0 * a2;
        if (d >= 0.0) {
            ix = 3;
            double sqd = std::sqrt(d);
            x[1] = 0.5 * (-a1 + sqd);
            x[2] = 0.5 * (-a1 - sqd);
        }
    }
    else {
        islv = 1;
        double Q = (3.0 * a2 - a1 * a1) / 9.0;
        double R = (9.0 * a1 * a2 - 27.0 * a3 - 2.0 * a1 * a1 * a1) / 54.0;
        double D = Q * Q * Q + R * R;

        if (D < -eps) {
            ix = 3;
            double arg = R / std::sqrt(-Q * Q * Q);
            arg = std::max(-1.0, std::min(arg, 1.0));
            double thet = expon * std::acos(arg);
            double coef = 2.0 * std::sqrt(-Q);
            x[0] = coef * std::cos(thet) - expon * a1;
            x[1] = coef * std::cos(thet + thet1 * pi) - expon * a1;
            x[2] = coef * std::cos(thet + thet2 * pi) - expon * a1;
        }
        else if (D <= eps) {
            ix = 2;
            double ssig = (R >= 0.0) ? 1.0 : -1.0;
            double s = ssig * std::pow(std::abs(R), expon);
            x[0] = 2.0 * s - expon * a1;
            x[1] = -s - expon * a1;
        }
        else {
            ix = 1;
            double sqd = std::sqrt(D);
            double ssig = (R + sqd >= 0.0) ? 1.0 : -1.0;
            double tsig = (R - sqd >= 0.0) ? 1.0 : -1.0;
            double s = ssig * std::pow(std::abs(R + sqd), expon);
            double t = tsig * std::pow(std::abs(R - sqd), expon);
            x[0] = s + t - expon * a1;
        }
    }

    root = 1.0e30;
    for (int i = 0; i < ix; ++i) {
        if (x[i] > 0.0) {
            root = std::min(root, x[i]);
            islv = 0;
        }
    }
}

void Solver::cal_claq(double cli, double hi, double& delt, State& state) {
    double a32 = state.xk31 * state.water / std::pow(state.gama[10], 2.0); // GAMA(11) -> index 10
    double om1 = cli;
    double om2 = hi;
    double bb = -(om1 + om2 + a32);
    double cc = om1 * om2;
    double dd = bb * bb - 4.0 * cc;
    if (dd >= 0.0) {
        double sqdd = std::sqrt(dd);
        double del1 = 0.5 * (-bb - sqdd);
        double del2 = 0.5 * (-bb + sqdd);
        if (del1 < 0.0 || del1 > hi || del1 > cli) {
            delt = 0.0;
        } else {
            delt = del1;
            return;
        }
        if (del2 < 0.0 || del2 > cli || del2 > hi) {
            delt = 0.0;
        } else {
            delt = del2;
        }
    } else {
        delt = 0.0;
    }
}

void Solver::cal_niaq(double no3i, double hi, double& delt, State& state) {
    double a42 = state.xk42 * state.water / std::pow(state.gama[9], 2.0); // GAMA(10) -> index 9
    double om1 = no3i;
    double om2 = hi;
    double bb = -(om1 + om2 + a42);
    double cc = om1 * om2;
    double dd = bb * bb - 4.0 * cc;
    if (dd >= 0.0) {
        double sqdd = std::sqrt(dd);
        double del1 = 0.5 * (-bb - sqdd);
        double del2 = 0.5 * (-bb + sqdd);
        if (del1 < 0.0 || del1 > hi || del1 > no3i) {
            delt = 0.0;
        } else {
            delt = del1;
            return;
        }
        if (del2 < 0.0 || del2 > no3i || del2 > hi) {
            delt = 0.0;
        } else {
            delt = del2;
        }
    } else {
        delt = 0.0;
    }
}

void Solver::isrp3f(const Input& input, State& state) {
    state.clear_errors();
    state.actmod = 3;

    // Adjust for too little ammonium and chloride
    state.w[2] = std::max(state.w[2], 1.0e-10); // NH3
    state.w[4] = std::max(state.w[4], 1.0e-10); // HCl

    if (state.w[0] + state.w[1] + state.w[3] <= 1.0e-10) {
        state.w[0] = 1.0e-10; // Na+
        state.w[1] = 1.0e-10; // SO4
    }

    double rest = 2.0 * state.w[1] + state.w[3] + state.w[4];
    if (state.w[0] > rest) {
        state.w[0] = (1.0 - 1e-6) * rest;
        state.push_error(50, "ISRP3F: SODIUM ADJUSTED");
    }

    double sulrat = (state.w[0] + state.w[2]) / state.w[1];
    double sodrat = state.w[0] / state.w[1];

    if (sulrat >= 2.0 && sodrat < 2.0) {
        state.scase = "G5";
        cal_cg5(input, state);
    }
    else if (sulrat >= 2.0 && sodrat >= 2.0) {
        state.scase = "H6";
        cal_ch6(input, state);
    }
    else if (sulrat >= 1.0 && sulrat < 2.0) {
        state.scase = "I6";
        cal_ci6(input, state);
        cal_cnha(input, state);
        cal_cnh3(input, state);
    }
    else {
        state.scase = "J3";
        cal_cj3(input, state);
        cal_cnha(input, state);
        cal_cnh3(input, state);
    }
}

//=======================================================================
// CASE G5 - SULFATE POOR, SODIUM POOR METASTABLE LIQUID Speciation
//=======================================================================
void Solver::cal_cg5(const Input& input, State& state) {
    state.calaou = true;
    double psi6lo = state.tiny;
    double psi6hi = state.w[4]; // Total Chloride

    double x1 = psi6lo;
    double y1 = funcg5a(x1, input, state);
    double eps = 1e-6;
    double x2 = x1;
    double y2 = y1;

    if (std::abs(y1) <= eps || state.w[4] <= state.tiny) {
        goto GOTO_50;
    }

    {
        int ndiv = 5;
        double dx = (psi6hi - psi6lo) / static_cast<double>(ndiv);
        bool sign_changed = false;

        for (int i = 1; i <= ndiv; ++i) {
            x2 = x1 + dx;
            state.rstgamp();
            y2 = funcg5a(x2, input, state);
            if ((y1 < 0.0 && y2 > 0.0) || (y1 > 0.0 && y2 < 0.0)) {
                sign_changed = true;
                break;
            }
            x1 = x2;
            y1 = y2;
        }

        if (!sign_changed) {
            if (std::abs(y2) > eps) {
                state.rstgamp();
                funcg5a(psi6lo, input, state);
            }
            goto GOTO_50;
        }
    }

    // Bisection
    {
        int maxit = 100;
        for (int i = 1; i <= maxit; ++i) {
            double x3 = 0.5 * (x1 + x2);
            state.rstgamp();
            double y3 = funcg5a(x3, input, state);

            if (((y1 < 0.0 && y3 <= 0.0) || (y1 > 0.0 && y3 >= 0.0)) == false) {
                y2 = y3;
                x2 = x3;
            } else {
                y1 = y3;
                x1 = x3;
            }

            if (std::abs(x2 - x1) <= eps * x1) {
                double final_x3 = 0.5 * (x1 + x2);
                state.rstgamp();
                funcg5a(final_x3, input, state);
                goto GOTO_50;
            }
        }
    }

    state.push_error(2, "CALCG5: NO CONVERGENCE");
    {
        double final_x3 = 0.5 * (x1 + x2);
        state.rstgamp();
        funcg5a(final_x3, input, state);
    }

GOTO_50:
    if (state.molal[1] > state.tiny && state.molal[5] > state.tiny) {
        double delta = 0.0;
        cal_chs4(state.molal[1], state.molal[5], 0.0, delta, state);
        state.molal[1] -= delta;
        state.molal[5] -= delta;
        state.molal[6] = delta;
    }
}

double Solver::funcg5a(double x, const Input& input, State& state) {
    state.frst   = true;
    state.calain = true;

    state.psi6 = x;
    state.psi1 = state.w[1]; // Total Sulfate
    state.psi2 = 0.0;
    state.psi3 = 0.0;
    state.psi7 = 0.5 * state.w[0];
    state.psi8 = 0.0;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double a4 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
        double a5 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
        double a6 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);

        if (state.w[3] > state.tiny) {
            state.psi5 = state.psi6 * state.w[3] / ((a6 / a5) * (state.w[4] - state.psi6) + state.psi6);
        } else {
            state.psi5 = state.tiny;
        }
        state.psi5 = std::max(state.psi5, state.tiny);

        state.psi4 = state.tiny;
        if (state.w[1] > state.tiny) {
            double bb = -(state.w[2] + state.psi6 + state.psi5 + 1.0 / a4);
            double cc = state.w[2] * (state.psi5 + state.psi6);
            double dd = std::max(bb * bb - 4.0 * cc, 0.0);
            state.psi4 = 0.5 * (-bb - std::sqrt(dd));
        }

        state.molal[0] = 2.0 * state.psi1;                     // Na+
        state.molal[2] = 2.0 * state.psi2 + state.psi4;        // NH4+
        state.molal[4] = state.psi6;                           // Cl-
        state.molal[5] = state.psi2 + state.psi1;              // SO4--
        state.molal[6] = 0.0;                                  // HSO4-
        state.molal[3] = state.psi5;                           // NO3-

        double smin = 2.0 * state.molal[5] + state.molal[3] + state.molal[4] - state.molal[0] - state.molal[2];
        double hi = 0.0;
        double ohi = 0.0;
        cal_cph(smin, hi, ohi, state);
        state.molal[1] = hi; // H+

        state.gnh3  = std::max(state.w[2] - state.psi4, state.tiny);
        state.ghno3 = std::max(state.w[3] - state.psi5, state.tiny);
        state.ghcl  = std::max(state.w[4] - state.psi6, state.tiny);

        state.cnacl   = 0.0;
        state.cnano3  = 0.0;
        state.cna2so4 = 0.0;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act3();
        } else {
            break;
        }
    }

    return state.molal[1] * state.molal[4] / state.ghcl / (state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0)) - 1.0;
}

//=======================================================================
// CASE H6 - SULFATE POOR, SODIUM RICH METASTABLE LIQUID Speciation
//=======================================================================
void Solver::cal_ch6(const Input& input, State& state) {
    state.calaou = true;
    double chi1 = state.w[1];
    double frna = std::max(state.w[0] - 2.0 * chi1, 0.0);
    double chi8 = std::min(frna, state.w[3]);
    double chi7 = std::min(std::max(frna - chi8, 0.0), state.w[4]);
    double chi6 = std::max(state.w[4] - chi7, 0.0);

    double psi6lo = state.tiny;
    double psi6hi = chi6 - state.tiny;

    double x1 = psi6lo;
    double y1 = funch6a(x1, input, state);
    double eps = 1e-6;
    double x2 = x1;
    double y2 = y1;

    if (std::abs(y1) <= eps || chi6 <= state.tiny) {
        goto GOTO_50;
    }

    {
        int ndiv = 5;
        double dx = (psi6hi - psi6lo) / static_cast<double>(ndiv);
        bool sign_changed = false;

        for (int i = 1; i <= ndiv; ++i) {
            x2 = x1 + dx;
            state.rstgamp();
            y2 = funch6a(x2, input, state);
            if ((y1 < 0.0 && y2 > 0.0) || (y1 > 0.0 && y2 < 0.0)) {
                sign_changed = true;
                break;
            }
            x1 = x2;
            y1 = y2;
        }

        if (!sign_changed) {
            if (std::abs(y2) > eps) {
                state.rstgamp();
                funch6a(psi6lo, input, state);
            }
            goto GOTO_50;
        }
    }

    // Bisection
    {
        int maxit = 100;
        for (int i = 1; i <= maxit; ++i) {
            double x3 = 0.5 * (x1 + x2);
            state.rstgamp();
            double y3 = funch6a(x3, input, state);

            if (((y1 < 0.0 && y3 <= 0.0) || (y1 > 0.0 && y3 >= 0.0)) == false) {
                y2 = y3;
                x2 = x3;
            } else {
                y1 = y3;
                x1 = x3;
            }

            if (std::abs(x2 - x1) <= eps * x1) {
                double final_x3 = 0.5 * (x1 + x2);
                state.rstgamp();
                funch6a(final_x3, input, state);
                goto GOTO_50;
            }
        }
    }

    state.push_error(2, "CALCH6: NO CONVERGENCE");
    {
        double final_x3 = 0.5 * (x1 + x2);
        state.rstgamp();
        funch6a(final_x3, input, state);
    }

GOTO_50:
    if (state.molal[1] > state.tiny && state.molal[5] > state.tiny) {
        double delta = 0.0;
        cal_chs4(state.molal[1], state.molal[5], 0.0, delta, state);
        state.molal[1] -= delta;
        state.molal[5] -= delta;
        state.molal[6] = delta;
    }
}

double Solver::funch6a(double x, const Input& input, State& state) {
    state.frst   = true;
    state.calain = true;

    state.psi6 = x;
    state.psi1 = state.w[1]; // Total Sulfate
    state.psi2 = 0.0;
    state.psi3 = 0.0;
    
    double frna = std::max(state.w[0] - 2.0 * state.psi1, 0.0);
    double chi8 = std::min(frna, state.w[3]);
    double chi5 = std::max(state.w[3] - chi8, 0.0);
    double chi7 = std::min(std::max(frna - chi8, 0.0), state.w[4]);
    double chi6 = std::max(state.w[4] - chi7, 0.0);

    state.psi7 = chi7;
    state.psi8 = chi8;
    state.psi4 = 0.0;
    state.psi5 = 0.0;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double a4 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
        double a5 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
        double a6 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);

        state.psi5 = chi5 * (state.psi6 + state.psi7) - (a6 / a5) * state.psi8 * (chi6 - state.psi6 - state.psi3);
        state.psi5 = state.psi5 / ((a6 / a5) * (chi6 - state.psi6 - state.psi3) + state.psi6 + state.psi7);
        state.psi5 = std::max(state.psi5, state.tiny);

        if (state.w[2] > state.tiny && state.water > state.tiny) {
            double bb = -(state.w[2] + state.psi6 + state.psi5 + 1.0 / a4);
            double cc = state.w[2] * (state.psi5 + state.psi6);
            double dd = bb * bb - 4.0 * cc;
            state.psi4 = 0.5 * (-bb - std::sqrt(dd));
            state.psi4 = std::min(state.psi4, state.w[2]);
        } else {
            state.psi4 = state.tiny;
        }

        state.molal[0] = state.psi8 + state.psi7 + 2.0 * state.psi1; // Na+
        state.molal[2] = state.psi4;                    // NH4+
        state.molal[4] = state.psi6 + state.psi7;              // Cl-
        state.molal[5] = state.psi2 + state.psi1;              // SO4--
        state.molal[6] = 0.0;                      // HSO4-
        state.molal[3] = state.psi5 + state.psi8;              // NO3-

        double smin = 2.0 * state.molal[5] + state.molal[3] + state.molal[4] - state.molal[0] - state.molal[2];
        double hi = 0.0;
        double ohi = 0.0;
        cal_cph(smin, hi, ohi, state);
        state.molal[1] = hi; // H+

        state.gnh3  = std::max(state.w[2] - state.psi4, state.tiny);
        state.ghno3 = std::max(chi5 - state.psi5, state.tiny);
        state.ghcl  = std::max(chi6 - state.psi6, state.tiny);

        state.cnacl   = std::max(chi7 - state.psi7, 0.0);
        state.cnano3  = std::max(chi8 - state.psi8, 0.0);
        state.cna2so4 = std::max(state.w[1] - state.psi1, 0.0);

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act3();
        } else {
            break;
        }
    }

    // Re-populate final speciation outside the loop (to prevent zero-out by cal_act3)
    state.molal[0] = state.psi8 + state.psi7 + 2.0 * state.psi1; // Na+
    state.molal[2] = state.psi4;                    // NH4+
    state.molal[4] = state.psi6 + state.psi7;              // Cl-
    state.molal[5] = state.psi2 + state.psi1;              // SO4--
    state.molal[6] = 0.0;                      // HSO4-
    state.molal[3] = state.psi5 + state.psi8;              // NO3-

    double smin = 2.0 * state.molal[5] + state.molal[3] + state.molal[4] - state.molal[0] - state.molal[2];
    double hi = 0.0, ohi = 0.0;
    cal_cph(smin, hi, ohi, state);
    state.molal[1] = hi;

    state.gnh3  = std::max(state.w[2] - state.psi4, state.tiny);
    state.ghno3 = std::max(chi5 - state.psi5, state.tiny);
    state.ghcl  = std::max(chi6 - state.psi6, state.tiny);

    return state.molal[2] * state.molal[4] / state.ghcl / state.gnh3 / (state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0)) / ( (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0) ) - 1.0;
}

//=======================================================================
// CASE I6 - SULFATE RICH, NO FREE ACID Speciation
//=======================================================================
void Solver::cal_ci6(const Input& input, State& state) {
    state.frst   = true;
    state.calain = true;
    state.calaou = true;

    // Dry composition material balance
    cal_ci1a(input, state);
    state.chi1 = state.cnh4hs4;
    state.chi2 = state.clc;
    state.chi3 = state.cnahso4;
    state.chi4 = state.cna2so4;
    state.chi5 = state.cnh42s4;

    state.psi1 = state.chi1;
    state.psi2 = state.chi2;
    state.psi3 = state.chi3;
    state.psi4 = state.chi4;
    state.psi5 = state.chi5;

    state.clc     = 0.0;
    state.cnh4hs4 = 0.0;
    state.cnahso4 = 0.0;
    state.cna2so4 = 0.0;
    state.cnh42s4 = 0.0;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double a6 = state.xk1 * state.water / state.gama[6] * std::pow(state.gama[7] / state.gama[6], 2.0);
        double bb = state.psi2 + state.psi4 + state.psi5 + a6;
        double cc = -a6 * (state.psi2 + state.psi3 + state.psi1);
        double dd = bb * bb - 4.0 * cc;
        state.psi6 = 0.5 * (-bb + std::sqrt(dd));

        state.molal[1] = state.psi6;                                       // H+
        state.molal[0] = 2.0 * state.psi4 + state.psi3;                    // Na+
        state.molal[2] = 3.0 * state.psi2 + 2.0 * state.psi5 + state.psi1; // NH4+
        state.molal[4] = 0.0;                                              // Cl-
        state.molal[5] = state.psi2 + state.psi4 + state.psi5 + state.psi6; // SO4--
        state.molal[6] = state.psi2 + state.psi3 + state.psi1 - state.psi6; // HSO4-
        state.molal[3] = 0.0;                                              // NO3-

        state.cal_cmr();

        if (!state.calain) break;
        state.cal_act3();
    }
}

void Solver::cal_ci1a(const Input& input, State& state) {
    state.cna2so4 = 0.5 * state.w[0];
    state.cnh4hs4 = 0.0;
    state.cnahso4 = 0.0;
    state.cnh42s4 = 0.0;
    double frso4 = std::max(state.w[1] - state.cna2so4, 0.0);

    if (state.w[2] <= frso4) {
        state.clc     = 0.0;
        state.cnh4hs4 = state.w[2];
        state.cnahso4 = frso4 - state.w[2];
        state.cna2so4 = state.cna2so4 - state.cnahso4;
    } else {
        double x = 2.0 * frso4 - state.w[2];
        double y = state.w[2] - frso4;
        if (x <= y) {
            state.clc     = x;
            state.cnh4hs4 = 0.0;
            state.cnh42s4 = y - x;
        } else {
            state.clc     = y;
            state.cnh4hs4 = x - y;
            state.cnh42s4 = 0.0;
        }
    }

    state.ghno3 = state.w[3];
    state.ghcl  = state.w[4];
    state.gnh3  = 0.0;
}

//=======================================================================
// CASE J3 - SULFATE RICH WITH FREE ACID Speciation
//=======================================================================
void Solver::cal_cj3(const Input& input, State& state) {
    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    double lamda = std::max(state.w[1] - state.w[2] - state.w[0], state.tiny);
    double chi1  = state.w[0];
    double chi2  = state.w[2];
    double psi1  = chi1;
    double psi2  = chi2;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        double a3 = state.xk1 * state.water / state.gama[6] * std::pow(state.gama[7] / state.gama[6], 2.0);
        double bb = a3 + lamda;
        double cc = -a3 * (lamda + psi1 + psi2);
        double dd = bb * bb - 4.0 * cc;
        double kapa = 0.5 * (-bb + std::sqrt(dd));

        state.molal[1] = lamda + kapa; // H+
        state.molal[0] = psi1;         // Na+
        state.molal[2] = psi2;         // NH4+
        state.molal[4] = 0.0;          // Cl-
        state.molal[5] = kapa;         // SO4--
        state.molal[6] = lamda + psi1 + psi2 - kapa; // HSO4-
        state.molal[3] = 0.0;          // NO3-

        state.cnahso4 = 0.0;
        state.cnh4hs4 = 0.0;

        state.cal_cmr();

        if (!state.calain) break;
        state.cal_act2();
    }
}

//=======================================================================
// RECOREGULATION: Dynamic pH / H+-OH Equilibria Solver (CALCPH Equivalent)
//=======================================================================
void Solver::cal_cph(double gg, double& hi, double& ohi, State& state) {
    double a7 = state.xkw * state.rh * state.water * state.water;
    double bb = gg;
    double cc = -a7;
    double dd = bb * bb - 4.0 * cc;
    hi = 0.5 * (bb + std::sqrt(dd));
    if (hi <= state.tiny) {
        double abb = std::abs(bb);
        double denm = (bb + abb) + 2.0 * a7 / abb;
        hi = 2.0 * a7 / denm;
    }
    ohi = a7 / hi;
}

void Solver::isrp4f(const Input& input, State& state) {
    state.clear_errors();
    state.scase = "4F";
    state.actmod = 4;

    // 1. Sequential adjustments for excess cations (replicates lines 265-330 of isofwd.f)
    double rest = 2.0 * state.w[1] + state.w[3] + state.w[4];
    if (state.w[0] + state.w[5] + state.w[6] + state.w[7] > rest) {
        double ccaso4_i  = std::min(state.w[1], state.w[5]);
        double frso4_i   = std::max(state.w[1] - ccaso4_i, 0.0);
        double ca_fr_i   = std::max(state.w[5] - ccaso4_i, 0.0);
        double ccano32_i = std::min(ca_fr_i, 0.5 * state.w[3]);
        ca_fr_i          = std::max(ca_fr_i - ccano32_i, 0.0);
        double no3_fr_i  = std::max(state.w[3] - 2.0 * ccano32_i, 0.0);
        double ccacl2_i  = std::min(ca_fr_i, 0.5 * state.w[4]);
        double cl_fr_i   = std::max(state.w[4] - 2.0 * ccacl2_i, 0.0);
        double rest1     = 2.0 * frso4_i + no3_fr_i + cl_fr_i;

        double cna2so4_i = std::min(frso4_i, 0.5 * state.w[0]);
        frso4_i          = std::max(frso4_i - cna2so4_i, 0.0);
        double na_fr_i   = std::max(state.w[0] - 2.0 * cna2so4_i, 0.0);
        double cnacl_i   = std::min(na_fr_i, cl_fr_i);
        na_fr_i          = std::max(na_fr_i - cnacl_i, 0.0);
        cl_fr_i          = std::max(cl_fr_i - cnacl_i, 0.0);
        double cnano3_i  = std::min(na_fr_i, no3_fr_i);
        double rest2     = 2.0 * frso4_i + no3_fr_i + cl_fr_i;

        double cmgso4_i  = std::min(frso4_i, state.w[7]);
        double fr_mg_i   = std::max(state.w[7] - cmgso4_i, 0.0);
        frso4_i          = std::max(frso4_i - cmgso4_i, 0.0);
        double cmgno32_i = std::min(fr_mg_i, 0.5 * no3_fr_i);
        fr_mg_i          = std::max(fr_mg_i - cmgno32_i, 0.0);
        no3_fr_i         = std::max(no3_fr_i - 2.0 * cmgno32_i, 0.0);
        double cmgcl2_i  = std::min(fr_mg_i, 0.5 * cl_fr_i);
        double rest3     = 2.0 * frso4_i + no3_fr_i + cl_fr_i;

        if (state.w[5] > rest) {
            state.w[5] = (1.0 - 1e-6) * rest;
            state.w[0] = 0.0;
            state.w[6] = 0.0;
            state.w[7] = 0.0;
            state.push_error(51, "ISRP4F: CA, NA, K, MG IN EXCESS");
        } else if (state.w[0] > rest1) {
            state.w[0] = (1.0 - 1e-6) * rest1;
            state.w[6] = 0.0;
            state.w[7] = 0.0;
            state.push_error(52, "ISRP4F: NA, K, MG IN EXCESS");
        } else if (state.w[7] > rest2) {
            state.w[7] = (1.0 - 1e-6) * rest2;
            state.w[6] = 0.0;
            state.push_error(53, "ISRP4F: K, MG IN EXCESS");
        } else if (state.w[6] > rest3) {
            state.w[6] = (1.0 - 1e-6) * rest3;
            state.push_error(54, "ISRP4F: K IN EXCESS");
        }
    }

    double so4rat  = (state.w[0] + state.w[2] + state.w[5] + state.w[6] + state.w[7]) / state.w[1];
    double crnarat = (state.w[0] + state.w[5] + state.w[6] + state.w[7]) / state.w[1];
    double crrat   = (state.w[5] + state.w[6] + state.w[7]) / state.w[1];

    if (so4rat >= 2.0 && crnarat < 2.0) {
        state.scase = "O7";
        cal_co7(input, state);
    }
    else if (so4rat >= 2.0 && crnarat >= 2.0 && crrat <= 2.0) {
        state.scase = "M8";
        cal_cm8(input, state);
    }
    else if (so4rat >= 2.0 && crnarat >= 2.0 && crrat > 2.0) {
        state.scase = "P13";
        cal_cp13(input, state);
    }
    else if (so4rat >= 1.0 && so4rat < 2.0) {
        state.scase = "L9";
        cal_cl9(input, state);
        cal_cnha(input, state);
        cal_cnh3(input, state);
    }
    else {
        state.scase = "K4";
        cal_ck4(input, state);
        cal_cnha(input, state);
        cal_cnh3(input, state);
    }
}

//=======================================================================
// CASE O7 - SULFATE POOR, DUST & SODIUM POOR METASTABLE LIQUID Speciation
//=======================================================================
void Solver::cal_co7(const Input& input, State& state) {
    state.calaou = true;
    state.chi9 = std::min(state.w[5], state.w[1]); // CCASO4
    double so4fr = std::max(state.w[1] - state.chi9, 0.0);
    state.chi7 = std::min(0.5 * state.w[6], so4fr); // CK2SO4
    so4fr        = std::max(so4fr - state.chi7, 0.0);
    state.chi1 = std::min(0.5 * state.w[0], so4fr); // NA2SO4
    so4fr        = std::max(so4fr - state.chi1, 0.0);
    state.chi8 = std::min(state.w[7], so4fr); // CMGSO4
    so4fr        = std::max(so4fr - state.chi8, 0.0);
    state.chi3 = 0.0;
    state.chi5 = state.w[3]; // Total Nitrate
    state.chi6 = state.w[4]; // Total Chloride
    state.chi2 = std::max(so4fr, 0.0);
    state.chi4 = std::max(state.w[2] - 2.0 * state.chi2, 0.0);

    state.psi1 = state.chi1;
    state.psi2 = state.chi2;
    state.psi3 = 0.0;
    state.psi4 = 0.0;
    state.psi5 = 0.0;
    state.psi6 = 0.0;
    state.psi7 = state.chi7;
    state.psi8 = state.chi8;

    double psi6lo = state.tiny;
    double psi6hi = state.chi6 - state.tiny;

    state.cal_cmr(); // Populate m0 lookups safely before dividing

    state.water = state.chi2 / state.m0[3] + state.chi1 / state.m0[1] + state.chi7 / state.m0[16] + state.chi8 / state.m0[20];
    state.water = std::max(state.water, state.tiny);

    // Initial values for bisection
    double x1 = psi6lo;
    double y1 = funco7(x1, input, state);

    if (state.chi6 <= state.tiny) {
        if (state.molal[1] > state.tiny && state.molal[5] > state.tiny) {
            double delta = 0.0;
            cal_chs4(state.molal[1], state.molal[5], 0.0, delta, state);
            state.molal[1] -= delta;
            state.molal[5] -= delta;
            state.molal[6] = delta;
        }
        return;
    }

    int ndiv = 5;
    double dx = (psi6hi - psi6lo) / static_cast<double>(ndiv);
    double x2 = x1;
    double y2 = y1;
    bool bracket_found = false;

    for (int i = 1; i <= ndiv; ++i) {
        x2 = x1 + dx;
        y2 = funco7(x2, input, state);
        if (y1 * y2 < 0.0) {
            bracket_found = true;
            break;
        }
        x1 = x2;
        y1 = y2;
    }

    if (!bracket_found) {
        if (std::abs(y2) > 1e-6) {
            state.rstgamp();
            y2 = funco7(psi6lo, input, state);
        }
        if (state.molal[1] > state.tiny && state.molal[5] > state.tiny) {
            double delta = 0.0;
            cal_chs4(state.molal[1], state.molal[5], 0.0, delta, state);
            state.molal[1] -= delta;
            state.molal[5] -= delta;
            state.molal[6] = delta;
        }
        return;
    }

    // Perform Bisection
    double eps = 1e-6;
    int maxit = 100;
    double x3 = 0.5 * (x1 + x2);
    double y3 = 0.0;
    bool converged = false;

    for (int it = 0; it < maxit; ++it) {
        x3 = 0.5 * (x1 + x2);
        state.rstgamp();
        y3 = funco7(x3, input, state);
        if (y1 * y3 <= 0.0) {
            y2 = y3;
            x2 = x3;
        } else {
            y1 = y3;
            x1 = x3;
        }
        if (std::abs(x2 - x1) <= eps * x1) {
            converged = true;
            break;
        }
    }

    if (!converged) {
        state.push_error(2, "CALCO7");
    }

    x3 = 0.5 * (x1 + x2);
    state.rstgamp();
    y3 = funco7(x3, input, state);

    if (state.molal[1] > state.tiny && state.molal[5] > state.tiny) {
        double delta = 0.0;
        cal_chs4(state.molal[1], state.molal[5], 0.0, delta, state);
        state.molal[1] -= delta;
        state.molal[5] -= delta;
        state.molal[6] = delta;
    }
}

double Solver::funco7(double x, const Input& input, State& state) {
    state.psi6 = x;
    state.frst = true;
    state.calain = true;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        state.a4 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
        state.a5 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
        state.a6 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);

        if (state.chi5 >= state.tiny) {
            state.psi5 = state.psi6 * state.chi5 / ((state.a6 / state.a5) * (state.chi6 - state.psi6) + state.psi6);
            state.psi5 = std::min(state.psi5, state.chi5);
        } else {
            state.psi5 = state.tiny;
        }

        if (state.w[2] > state.tiny) {
            double bb = -(state.chi4 + state.psi6 + state.psi5 + 1.0 / state.a4);
            double cc = state.chi4 * (state.psi5 + state.psi6) - 2.0 * state.psi2 / state.a4;
            double dd = std::max(bb * bb - 4.0 * cc, 0.0);
            state.psi4 = 0.5 * (-bb - std::sqrt(dd));
            state.psi4 = std::max(std::min(state.psi4, state.chi4), 0.0);
        } else {
            state.psi4 = state.tiny;
        }

        state.molal[0] = 2.0 * state.psi1;                     // Na+
        state.molal[2] = 2.0 * state.psi2 + state.psi4;        // NH4+
        state.molal[4] = state.psi6;                           // Cl-
        state.molal[5] = state.psi1 + state.psi2 + state.psi7 + state.psi8; // SO4--
        state.molal[6] = 0.0;                                  // HSO4-
        state.molal[3] = state.psi5;                           // NO3-
        state.molal[7] = 0.0;                                  // Ca++
        state.molal[8] = 2.0 * state.psi7;                     // K+
        state.molal[9] = state.psi8;                           // Mg++

        double smin = 2.0 * state.molal[5] + state.molal[3] + state.molal[4] 
                    - state.molal[0] - state.molal[2] - state.molal[8] - 2.0 * state.molal[9];
        double hi = 0.0, ohi = 0.0;
        cal_cph(smin, hi, ohi, state);
        state.molal[1] = hi;

        state.gnh3  = std::max(state.chi4 - state.psi4, state.tiny);
        state.ghno3 = std::max(state.chi5 - state.psi5, state.tiny);
        state.ghcl  = std::max(state.chi6 - state.psi6, state.tiny);

        state.cna2so4 = 0.0;
        state.cnh42s4 = 0.0;
        state.cnh4no3 = 0.0;
        state.cnh4cl  = 0.0;
        state.ck2so4  = 0.0;
        state.cmgso4  = 0.0;
        state.ccaso4  = state.chi9;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act4();
        } else {
            break;
        }
    }

    return state.molal[1] * state.molal[4] / state.ghcl / state.a6 - 1.0;
}

//=======================================================================
// CASE M8: Sulfate poor, dust poor Metastable liquid speciation
//=======================================================================
void Solver::cal_cm8(const Input& input, State& state) {
    double eps = 1e-6;
    int ndiv = 5;
    int maxit = 100;

    state.calaou = true;
    state.chi11  = std::min(state.w[5], state.w[1]);             // CCASO4
    double so4fr  = std::max(state.w[1] - state.chi11, 0.0);
    double cafr   = std::max(state.w[5] - state.chi11, 0.0);
    state.chi9   = std::min(0.5 * state.w[6], so4fr);             // CK2SO4
    double frk    = std::max(state.w[6] - 2.0 * state.chi9, 0.0);
    so4fr  = std::max(so4fr - state.chi9, 0.0);
    state.chi10  = std::min(state.w[7], so4fr);                  // CMGSO4
    double frmg   = std::max(state.w[7] - state.chi10, 0.0);
    so4fr  = std::max(so4fr - state.chi10, 0.0);
    state.chi1   = std::max(so4fr, 0.0);                         // CNA2SO4
    state.chi2   = 0.0;                                          // CNH42S4
    state.chi3   = 0.0;                                          // CNH4CL
    double frna   = std::max(state.w[0] - 2.0 * state.chi1, 0.0);
    state.chi8   = std::min(frna, state.w[3]);                   // CNANO3
    state.chi4   = state.w[2];                                   // NH3(g)
    state.chi5   = std::max(state.w[3] - state.chi8, 0.0);       // HNO3(g)
    state.chi7   = std::min(std::max(frna - state.chi8, 0.0), state.w[4]); // CNACL
    state.chi6   = std::max(state.w[4] - state.chi7, 0.0);       // HCL(g)

    state.cal_cmr(); // Populate m0 lookups safely before dividing

    state.water = state.chi1 / state.m0[1] + state.chi7 / state.m0[0] + state.chi8 / state.m0[2] +
                  state.chi9 / state.m0[16] + state.chi10 / state.m0[20];
    if (std::isnan(state.water)) {
        std::printf("NAN FOUND inside cal_cm8: chi1=%e, m0[1]=%e, chi7=%e, m0[0]=%e, chi8=%e, m0[2]=%e, chi9=%e, m0[16]=%e, chi10=%e, m0[20]=%e\n",
                    state.chi1, state.m0[1], state.chi7, state.m0[0], state.chi8, state.m0[2], state.chi9, state.m0[16], state.chi10, state.m0[20]);
    }
    state.water = std::max(state.water, state.tiny);

    double psi6lo = state.tiny;
    double psi6hi = state.chi6 - state.tiny;

    double x1 = psi6lo;
    double y1 = funcm8(x1, input, state);

    if (std::abs(y1) <= eps || state.chi6 <= state.tiny) {
        // Early convergence or negligible HCl
    } else {
        // ROOT TRACKING
        double dx = (psi6hi - psi6lo) / static_cast<double>(ndiv);
        double x2 = x1;
        double y2 = y1;
        bool bracket_found = false;

        for (int i = 1; i <= ndiv; ++i) {
            x2 = x1 + dx;
            y2 = funcm8(x2, input, state);
            if (y1 * y2 < 0.0) {
                bracket_found = true;
                break;
            }
            x1 = x2;
            y1 = y2;
        }

        if (!bracket_found) {
            if (std::abs(y2) > eps) {
                state.rstgamp();
                y2 = funcm8(psi6lo, input, state);
            }
        } else {
            // PERFORM BISECTION
            double x3 = 0.5 * (x1 + x2);
            bool converged = false;
            for (int it = 0; it < maxit; ++it) {
                x3 = 0.5 * (x1 + x2);
                state.rstgamp();
                double y3 = funcm8(x3, input, state);
                if (y1 * y3 <= 0.0) {
                    y2 = y3;
                    x2 = x3;
                } else {
                    y1 = y3;
                    x1 = x3;
                }
                if (std::abs(x2 - x1) <= eps * x1) {
                    converged = true;
                    break;
                }
            }

            if (!converged) {
                state.push_error(2, "CALCM8");
            }

            x3 = 0.5 * (x1 + x2);
            state.rstgamp();
            funcm8(x3, input, state);
        }
    }

    if (state.molal[1] > state.tiny && state.molal[5] > state.tiny) {
        double delta = 0.0;
        cal_chs4(state.molal[1], state.molal[5], 0.0, delta, state);
        state.molal[1] -= delta;
        state.molal[5] -= delta;
        state.molal[6] = delta;
    }
}

double Solver::funcm8(double x, const Input& input, State& state) {
    state.psi6   = x;
    state.psi1   = state.chi1;
    state.psi2   = 0.0;
    state.psi3   = 0.0;
    state.psi7   = state.chi7;
    state.psi8   = state.chi8;
    state.psi9   = state.chi9;
    state.psi10  = state.chi10;
    state.psi11  = 0.0;
    state.frst   = true;
    state.calain = true;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        state.a4 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
        state.a5 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
        state.a6 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);

        // Calculate dissociation quantities
        state.psi5 = state.chi5 * (state.psi6 + state.psi7) - (state.a6 / state.a5) * state.psi8 * (state.chi6 - state.psi6 - state.psi3);
        state.psi5 = state.psi5 / ((state.a6 / state.a5) * (state.chi6 - state.psi6 - state.psi3) + state.psi6 + state.psi7);
        state.psi5 = std::min(std::max(state.psi5, state.tiny), state.chi5);

        if (state.w[2] > state.tiny && state.water > state.tiny) {
            double bb = -(state.chi4 + state.psi6 + state.psi5 + 1.0 / state.a4);
            double cc = state.chi4 * (state.psi5 + state.psi6);
            double dd = std::max(bb * bb - 4.0 * cc, 0.0);
            state.psi4 = 0.5 * (-bb - std::sqrt(dd));
            state.psi4 = std::min(std::max(state.psi4, 0.0), state.chi4);
        } else {
            state.psi4 = state.tiny;
        }

        // Calculate speciation (0-based indexing mapping 1-to-1 to Fortran MOLAL)
        state.molal[0] = state.psi8 + state.psi7 + 2.0 * state.psi1;                     // Na+
        state.molal[2] = state.psi4;                                                     // NH4+
        state.molal[4] = state.psi6 + state.psi7;                                         // Cl-
        state.molal[5] = state.psi2 + state.psi1 + state.psi9 + state.psi10;              // SO4--
        state.molal[6] = 0.0;                                                            // HSO4-
        state.molal[3] = state.psi5 + state.psi8;                                         // NO3-
        state.molal[7] = state.psi11;                                                    // Ca++
        state.molal[8] = 2.0 * state.psi9;                                               // K+
        state.molal[9] = state.psi10;                                                    // Mg++

        double smin = 2.0 * state.molal[5] + state.molal[3] + state.molal[4] - state.molal[0] - state.molal[2]
                    - state.molal[8] - 2.0 * state.molal[9];
        double hi = 0.0, ohi = 0.0;
        cal_cph(smin, hi, ohi, state);

        state.molal[1] = hi; // H+

        state.gnh3  = std::max(state.chi4 - state.psi4, state.tiny);
        state.ghno3 = std::max(state.chi5 - state.psi5, state.tiny);
        state.ghcl  = std::max(state.chi6 - state.psi6, state.tiny);

        state.cna2so4 = 0.0;
        state.cnh42s4 = 0.0;
        state.cnh4no3 = 0.0;
        state.cnacl   = 0.0;
        state.cnano3  = 0.0;
        state.cna2so4 = 0.0;
        state.ck2so4  = 0.0;
        state.cmgso4  = 0.0;
        state.ccaso4  = state.chi11;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act4();
        } else {
            break;
        }
    }

    // Re-populate final speciation outside the loop (to prevent zero-out by cal_act4)
    state.molal[0] = state.psi8 + state.psi7 + 2.0 * state.psi1;                     // Na+
    state.molal[2] = state.psi4;                                                     // NH4+
    state.molal[4] = state.psi6 + state.psi7;                                         // Cl-
    state.molal[5] = state.psi2 + state.psi1 + state.psi9 + state.psi10;              // SO4--
    state.molal[6] = 0.0;                                                            // HSO4-
    state.molal[3] = state.psi5 + state.psi8;                                         // NO3-
    state.molal[7] = state.psi11;                                                    // Ca++
    state.molal[8] = 2.0 * state.psi9;                                               // K+
    state.molal[9] = state.psi10;                                                    // Mg++

    double smin = 2.0 * state.molal[5] + state.molal[3] + state.molal[4] - state.molal[0] - state.molal[2]
                - state.molal[8] - 2.0 * state.molal[9];
    double hi = 0.0, ohi = 0.0;
    cal_cph(smin, hi, ohi, state);
    state.molal[1] = hi;

    state.gnh3  = std::max(state.chi4 - state.psi4, state.tiny);
    state.ghno3 = std::max(state.chi5 - state.psi5, state.tiny);
    state.ghcl  = std::max(state.chi6 - state.psi6, state.tiny);

    return state.molal[1] * state.molal[4] / state.ghcl / state.a6 - 1.0;
}

//=======================================================================
// CASE P13: Sulfate poor, dust rich Metastable liquid speciation
//=======================================================================
void Solver::cal_cp13(const Input& input, State& state) {
    double eps = 1e-6;
    int ndiv = 5;
    int maxit = 100;

    state.calaou  = true;
    state.chi11   = std::min(state.w[1], state.w[5]);                    // CCASO4
    double frca    = std::max(state.w[5] - state.chi11, 0.0);
    double frso4   = std::max(state.w[1] - state.chi11, 0.0);
    state.chi9    = std::min(frso4, 0.5 * state.w[6]);                    // CK2SO4
    double frk     = std::max(state.w[6] - 2.0 * state.chi9, 0.0);
    frso4   = std::max(frso4 - state.chi9, 0.0);
    state.chi10   = frso4;                                                // CMGSO4
    double frmg    = std::max(state.w[7] - state.chi10, 0.0);
    state.chi7    = std::min(state.w[0], state.w[4]);                     // CNACL
    double frna    = std::max(state.w[0] - state.chi7, 0.0);
    double frcl    = std::max(state.w[4] - state.chi7, 0.0);
    state.chi12   = std::min(frca, 0.5 * state.w[3]);                     // CCANO32
    frca    = std::max(frca - state.chi12, 0.0);
    double frno3   = std::max(state.w[3] - 2.0 * state.chi12, 0.0);
    state.chi17   = std::min(frca, 0.5 * frcl);                           // CCACL2
    frca    = std::max(frca - state.chi17, 0.0);
    frcl    = std::max(frcl - 2.0 * state.chi17, 0.0);
    state.chi15   = std::min(frmg, 0.5 * frno3);                          // CMGNO32
    frmg    = std::max(frmg - state.chi15, 0.0);
    frno3   = std::max(frno3 - 2.0 * state.chi15, 0.0);
    state.chi16   = std::min(frmg, 0.5 * frcl);                           // CMGCL2
    frmg    = std::max(frmg - state.chi16, 0.0);
    frcl    = std::max(frcl - 2.0 * state.chi16, 0.0);
    state.chi8    = std::min(frna, frno3);                                // CNANO3
    frna    = std::max(frna - state.chi8, 0.0);
    frno3   = std::max(frno3 - state.chi8, 0.0);
    state.chi14   = std::min(frk, frcl);                                  // CKCL
    frk     = std::max(frk - state.chi14, 0.0);
    frcl    = std::max(frcl - state.chi14, 0.0);
    state.chi13   = std::min(frk, frno3);                                 // CKNO3
    frk     = std::max(frk - state.chi13, 0.0);
    frno3   = std::max(frno3 - state.chi13, 0.0);

    state.chi5    = frno3;                                                // HNO3(g)
    state.chi6    = frcl;                                                 // HCL(g)
    state.chi4    = state.w[2];                                           // NH3(g)

    state.chi3    = 0.0;                                                  // CNH4CL
    state.chi1    = 0.0;
    state.chi2    = 0.0;

    state.cal_cmr(); // Populate m0 lookups safely before dividing

    state.water = state.chi9 / state.m0[16] + state.chi10 / state.m0[20] + state.chi7 / state.m0[0] +
                  state.chi12 / state.m0[14] + state.chi17 / state.m0[15] + state.chi15 / state.m0[21] +
                  state.chi16 / state.m0[22] + state.chi8 / state.m0[2] + state.chi14 / state.m0[19] +
                  state.chi13 / state.m0[18];
    state.water = std::max(state.water, state.tiny);

    double psi6lo = state.tiny;
    double psi6hi = state.chi6 - state.tiny;

    double x1 = psi6lo;
    double y1 = funcp13(x1, input, state);

    if (std::abs(y1) <= eps || state.chi6 <= state.tiny) {
        // Early convergence or negligible HCl
    } else {
        // ROOT TRACKING
        double dx = (psi6hi - psi6lo) / static_cast<double>(ndiv);
        double x2 = x1;
        double y2 = y1;
        bool bracket_found = false;

        for (int i = 1; i <= ndiv; ++i) {
            x2 = x1 + dx;
            y2 = funcp13(x2, input, state);
            if (y1 * y2 < 0.0) {
                bracket_found = true;
                break;
            }
            x1 = x2;
            y1 = y2;
        }

        if (!bracket_found) {
            if (std::abs(y2) > eps) {
                state.rstgamp();
                y2 = funcp13(psi6lo, input, state);
            }
        } else {
            // PERFORM BISECTION
            double x3 = 0.5 * (x1 + x2);
            bool converged = false;
            for (int it = 0; it < maxit; ++it) {
                x3 = 0.5 * (x1 + x2);
                state.rstgamp();
                double y3 = funcp13(x3, input, state);
                if (y1 * y3 <= 0.0) {
                    y2 = y3;
                    x2 = x3;
                } else {
                    y1 = y3;
                    x1 = x3;
                }
                if (std::abs(x2 - x1) <= eps * x1) {
                    converged = true;
                    break;
                }
            }

            if (!converged) {
                state.push_error(2, "CALCP13");
            }

            x3 = 0.5 * (x1 + x2);
            state.rstgamp();
            funcp13(x3, input, state);
        }
    }

    if (state.molal[1] > state.tiny && state.molal[5] > state.tiny) {
        double delta = 0.0;
        cal_chs4(state.molal[1], state.molal[5], 0.0, delta, state);
        state.molal[1] -= delta;
        state.molal[5] -= delta;
        state.molal[6] = delta;
    }
}

double Solver::funcp13(double x, const Input& input, State& state) {
    state.psi6   = x;
    state.psi1   = 0.0;
    state.psi2   = 0.0;
    state.psi3   = 0.0;
    state.psi4   = 0.0;
    state.psi7   = state.chi7;
    state.psi8   = state.chi8;
    state.psi9   = state.chi9;
    state.psi10  = state.chi10;
    state.psi11  = 0.0;
    state.psi12  = state.chi12;
    state.psi13  = state.chi13;
    state.psi14  = state.chi14;
    state.psi15  = state.chi15;
    state.psi16  = state.chi16;
    state.psi17  = state.chi17;
    state.frst   = true;
    state.calain = true;

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        state.a4 = (state.xk2 / state.xkw) * state.r * state.temp * std::pow(state.gama[9] / state.gama[4], 2.0);
        state.a5 = state.xk4 * state.r * state.temp * std::pow(state.water / state.gama[9], 2.0);
        state.a6 = state.xk3 * state.r * state.temp * std::pow(state.water / state.gama[10], 2.0);

        // Calculate dissociation quantities
        state.psi5 = state.chi5 * (state.psi6 + state.psi7 + state.psi14 + 2.0 * state.psi16 + 2.0 * state.psi17) -
                     (state.a6 / state.a5) * (state.psi8 + 2.0 * state.psi12 + state.psi13 + 2.0 * state.psi15) * (state.chi6 - state.psi6);
        state.psi5 = state.psi5 / ((state.a6 / state.a5) * (state.chi6 - state.psi6) + state.psi6 + state.psi7 + state.psi14 +
                     2.0 * state.psi16 + 2.0 * state.psi17);
        state.psi5 = std::min(std::max(state.psi5, state.tiny), state.chi5);

        if (state.w[2] > state.tiny && state.water > state.tiny) {
            double bb = -(state.chi4 + state.psi6 + state.psi5 + 1.0 / state.a4);
            double cc = state.chi4 * (state.psi5 + state.psi6);
            double dd = std::max(bb * bb - 4.0 * cc, 0.0);
            state.psi4 = 0.5 * (-bb - std::sqrt(dd));
            state.psi4 = std::min(std::max(state.psi4, 0.0), state.chi4);
        } else {
            state.psi4 = state.tiny;
        }

        // Calculate speciation
        state.molal[0] = state.psi8 + state.psi7;                                                         // Na+
        state.molal[2] = state.psi4;                                                                      // NH4+
        state.molal[4] = state.psi6 + state.psi7 + state.psi14 + 2.0 * state.psi16 + 2.0 * state.psi17;   // Cl-
        state.molal[5] = state.psi9 + state.psi10;                                                        // SO4--
        state.molal[6] = 0.0;                                                                             // HSO4-
        state.molal[3] = state.psi5 + state.psi8 + 2.0 * state.psi12 + state.psi13 + 2.0 * state.psi15;   // NO3-
        state.molal[7] = state.psi11 + state.psi12 + state.psi17;                                         // Ca++
        state.molal[8] = 2.0 * state.psi9 + state.psi13 + state.psi14;                                    // K+
        state.molal[9] = state.psi10 + state.psi15 + state.psi16;                                         // Mg++

        double smin = 2.0 * state.molal[5] + state.molal[3] + state.molal[4] - state.molal[0] - state.molal[2]
                    - state.molal[8] - 2.0 * state.molal[9] - 2.0 * state.molal[7];
        double hi = 0.0, ohi = 0.0;
        cal_cph(smin, hi, ohi, state);
        state.molal[1] = hi;

        state.gnh3  = std::max(state.chi4 - state.psi4, state.tiny);
        state.ghno3 = std::max(state.chi5 - state.psi5, state.tiny);
        state.ghcl  = std::max(state.chi6 - state.psi6, state.tiny);

        state.cnh4no3 = 0.0;
        state.cnh4cl  = 0.0;
        state.cnacl   = 0.0;
        state.cnano3  = 0.0;
        state.ck2so4  = 0.0;
        state.cmgso4  = 0.0;
        state.ccaso4  = state.chi11;
        state.ccano32  = 0.0;
        state.ckno3   = 0.0;
        state.ckcl    = 0.0;
        state.cmgno32 = 0.0;
        state.cmgcl2  = 0.0;
        state.ccacl2  = 0.0;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act4();
        } else {
            break;
        }
    }

    // Re-populate final speciation outside the loop (to prevent zero-out by cal_act4)
    state.molal[0] = state.psi8 + state.psi7;                                                         // Na+
    state.molal[2] = state.psi4;                                                                      // NH4+
    state.molal[4] = state.psi6 + state.psi7 + state.psi14 + 2.0 * state.psi16 + 2.0 * state.psi17;   // Cl-
    state.molal[5] = state.psi9 + state.psi10;                                                        // SO4--
    state.molal[6] = 0.0;                                                                             // HSO4-
    state.molal[3] = state.psi5 + state.psi8 + 2.0 * state.psi12 + state.psi13 + 2.0 * state.psi15;   // NO3-
    state.molal[7] = state.psi11 + state.psi12 + state.psi17;                                         // Ca++
    state.molal[8] = 2.0 * state.psi9 + state.psi13 + state.psi14;                                    // K+
    state.molal[9] = state.psi10 + state.psi15 + state.psi16;                                         // Mg++

    double smin = 2.0 * state.molal[5] + state.molal[3] + state.molal[4] - state.molal[0] - state.molal[2]
                - state.molal[8] - 2.0 * state.molal[9] - 2.0 * state.molal[7];
    double hi = 0.0, ohi = 0.0;
    cal_cph(smin, hi, ohi, state);
    state.molal[1] = hi;

    state.gnh3  = std::max(state.chi4 - state.psi4, state.tiny);
    state.ghno3 = std::max(state.chi5 - state.psi5, state.tiny);
    state.ghcl  = std::max(state.chi6 - state.psi6, state.tiny);

    return state.molal[1] * state.molal[4] / state.ghcl / state.a6 - 1.0;
}

//=======================================================================
// CASE L9: Sulfate rich, no free acid (1.0 <= SULRAT < 2.0) Speciation
//=======================================================================
void Solver::cal_cl9(const Input& input, State& state) {
    // 1. Find dry composition
    cal_cl1a(input, state);

    // 2. Setup parameters (Common SOLUT equivalents)
    state.chi1 = state.cnh4hs4;
    state.chi2 = state.clc;
    state.chi3 = state.cnahso4;
    state.chi4 = state.cna2so4;
    state.chi5 = state.cnh42s4;
    state.chi6 = state.ck2so4;
    state.chi7 = state.cmgso4;
    state.chi8 = state.ckhso4;

    state.psi1 = state.chi1;
    state.psi2 = state.chi2;
    state.psi3 = state.chi3;
    state.psi4 = state.chi4;
    state.psi5 = state.chi5;
    state.psi6 = state.chi6;
    state.psi7 = state.chi7;
    state.psi8 = state.chi8;

    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    // 3. Solve equations with iterations for Activity Coefficients
    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        state.a9 = state.xk1 * state.water / state.gama[6] * std::pow(state.gama[7] / state.gama[6], 2.0);

        // Calculate dissociation quantities
        double bb = state.psi7 + state.psi6 + state.psi5 + state.psi4 + state.psi2 + state.a9; // LAMDA
        double cc = -state.a9 * (state.psi8 + state.psi1 + state.psi2 + state.psi3);
        double dd = std::max(bb * bb - 4.0 * cc, 0.0);
        double lamda = 0.5 * (-bb + std::sqrt(dd));
        lamda = std::min(std::max(lamda, state.tiny), state.psi8 + state.psi3 + state.psi2 + state.psi1);

        // Save concentrations in molal array (0-based indexing)
        state.molal[0] = 2.0 * state.psi4 + state.psi3;                                  // Na+
        state.molal[1] = lamda;                                                         // H+
        state.molal[2] = 3.0 * state.psi2 + 2.0 * state.psi5 + state.psi1;              // NH4+
        state.molal[4] = 0.0;                                                           // Cl-
        state.molal[5] = state.psi2 + state.psi4 + state.psi5 + state.psi6 + state.psi7 + lamda; // SO4--
        state.molal[6] = state.psi2 + state.psi3 + state.psi1 + state.psi8 - lamda;         // HSO4-
        state.molal[3] = 0.0;                                                           // NO3-
        state.molal[7] = 0.0;                                                           // Ca++
        state.molal[8] = state.psi8 + 2.0 * state.psi6;                                 // K+
        state.molal[9] = state.psi7;                                                    // Mg++

        state.clc     = 0.0;
        state.cnahso4 = 0.0;
        state.cna2so4 = 0.0;
        state.cnh42s4 = 0.0;
        state.cnh4hs4 = 0.0;
        state.ck2so4  = 0.0;
        state.cmgso4  = 0.0;
        state.ckhso4  = 0.0;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act4();
        } else {
            break;
        }
    }
}

void Solver::cal_cl1a(const Input& input, State& state) {
    state.ccaso4  = std::min(state.w[5], state.w[1]);                           // CCASO4
    double frso4   = std::max(state.w[1] - state.ccaso4, 0.0);
    double cafr    = std::max(state.w[5] - state.ccaso4, 0.0);
    state.ck2so4  = std::min(0.5 * state.w[6], frso4);                          // CK2SO4
    double frk     = std::max(state.w[6] - 2.0 * state.ck2so4, 0.0);
    frso4          = std::max(frso4 - state.ck2so4, 0.0);
    state.cna2so4 = std::min(0.5 * state.w[0], frso4);                          // CNA2SO4
    double frna    = std::max(state.w[0] - 2.0 * state.cna2so4, 0.0);
    frso4          = std::max(frso4 - state.cna2so4, 0.0);
    state.cmgso4  = std::min(state.w[7], frso4);                                // CMGSO4
    double frmg    = std::max(state.w[7] - state.cmgso4, 0.0);
    frso4          = std::max(frso4 - state.cmgso4, 0.0);

    state.cnh4hs4 = 0.0;
    state.cnahso4 = 0.0;
    state.cnh42s4 = 0.0;
    state.ckhso4  = 0.0;

    state.clc     = std::min(state.w[2] / 3.0, frso4 / 2.0);
    frso4         = std::max(frso4 - 2.0 * state.clc, 0.0);
    double frnh4  = std::max(state.w[2] - 3.0 * state.clc, 0.0);

    if (frso4 <= state.tiny) {
        state.clc     = std::max(state.clc - frnh4, 0.0);
        state.cnh42s4 = 2.0 * frnh4;
    } else if (frnh4 <= state.tiny) {
        state.cnh4hs4 = 3.0 * std::min(frso4, state.clc);
        state.clc     = std::max(state.clc - frso4, 0.0);

        if (state.cna2so4 > state.tiny) {
            frso4          = std::max(frso4 - state.cnh4hs4 / 3.0, 0.0);
            state.cnahso4 = 2.0 * frso4;
            state.cna2so4  = std::max(state.cna2so4 - frso4, 0.0);
        }
        if (state.ck2so4 > state.tiny) {
            frso4         = std::max(frso4 - state.cnh4hs4 / 3.0, 0.0);
            state.ckhso4  = 2.0 * frso4;
            state.ck2so4  = std::max(state.ck2so4 - frso4, 0.0);
        }
    }

    state.ghno3 = state.w[3];
    state.ghcl  = state.w[4];
    state.gnh3  = 0.0;
}

//=======================================================================
// CASE K4: Sulfate super rich, free acid (SO4RAT < 1.0) Speciation
//=======================================================================
void Solver::cal_ck4(const Input& input, State& state) {
    state.calaou = true;
    state.frst   = true;
    state.calain = true;

    state.chi1 = state.w[2]; // Total NH4 initially as NH4HSO4
    state.chi2 = state.w[0]; // Total NA initially as NaHSO4
    state.chi3 = state.w[6]; // Total K initially as KHSO4
    state.chi4 = state.w[7]; // Total Mg initially as MgSO4

    double lamda = std::max(state.w[1] - state.w[2] - state.w[0] - state.w[5] - state.w[6] - state.w[7], state.tiny); // FREE H2SO4
    state.psi1 = state.chi1; // ALL NH4HSO4 DELIQUESCED
    state.psi2 = state.chi2; // ALL NaHSO4 DELIQUESCED
    state.psi3 = state.chi3; // ALL KHSO4 DELIQUESCED
    state.psi4 = state.chi4; // ALL MgSO4 DELIQUESCED

    int nsweep = 4;
    for (int sweep = 0; sweep < nsweep; ++sweep) {
        state.a4 = state.xk1 * state.water / state.gama[6] * std::pow(state.gama[7] / state.gama[6], 2.0);

        double bb   = state.a4 + lamda + state.psi4; // KAPA
        double cc   = -state.a4 * (lamda + state.psi3 + state.psi2 + state.psi1) + lamda * state.psi4;
        double dd   = std::max(bb * bb - 4.0 * cc, 0.0);
        double kapa = 0.5 * (-bb + std::sqrt(dd));

        // Save concentrations in molal array (0-based indexing)
        state.molal[0] = state.psi2;                                                     // Na+
        state.molal[1] = std::max(lamda + kapa, state.tiny);                             // H+
        state.molal[2] = state.psi1;                                                     // NH4+
        state.molal[4] = 0.0;                                                           // Cl-
        state.molal[5] = std::max(kapa + state.psi4, 0.0);                               // SO4--
        state.molal[6] = std::max(lamda + state.psi1 + state.psi2 + state.psi3 - kapa, 0.0); // HSO4-
        state.molal[3] = 0.0;                                                           // NO3-
        state.molal[7] = 0.0;                                                           // Ca++
        state.molal[8] = state.psi3;                                                    // K+
        state.molal[9] = state.psi4;                                                    // Mg++

        state.cnh4hs4 = 0.0;
        state.cnahso4 = 0.0;
        state.ckhso4  = 0.0;
        state.ccaso4  = state.w[5];
        state.cmgso4  = 0.0;

        state.cal_cmr();

        if (state.frst && state.calaou || !state.frst && state.calain) {
            state.cal_act4();
        } else {
            break;
        }
    }
}

} // namespace Isorropia
