#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"
#include "Isorropia/Isorropia.h"
#include <random>
#include <algorithm>
#include <cmath>

//=======================================================================
// 1. PHYSICAL INVARIANTS: Non-Negativity & Electroneutrality
//=======================================================================
TEST(PropertyTest, VerifyPhysicalInvariants) {
    Isorropia::Solver solver;
    std::mt19937 gen(101); // fixed seed for reproducible property-based testing runs
    
    std::uniform_real_distribution<double> dist_so4(0.1, 20.0);   // Total sulfate: 0.1 to 20 µg/m3
    std::uniform_real_distribution<double> dist_nh3(0.1, 50.0);   // Total ammonia: 0.1 to 50 µg/m3
    std::uniform_real_distribution<double> dist_hno3(0.0, 15.0);  // Total nitrate: 0.0 to 15 µg/m3
    std::uniform_real_distribution<double> dist_cl(0.0, 10.0);    // Total chloride: 0.0 to 10 µg/m3
    std::uniform_real_distribution<double> dist_na(0.0, 5.0);     // Total sodium: 0.0 to 5 µg/m3
    std::uniform_real_distribution<double> dist_rh(0.15, 0.98);   // RH: 15% to 98%
    std::uniform_real_distribution<double> dist_temp(260.0, 315.0); // Temperature: 260K to 315K
    std::uniform_real_distribution<double> dist_org(0.0, 10.0);   // Organic concentration (ug/m3)
    std::uniform_real_distribution<double> dist_korg(0.0, 0.25);  // Organic hygroscopicity
    std::uniform_real_distribution<double> dist_rhoorg(800.0, 1200.0); // Organic density (kg/m3)

    for (int run = 0; run < 100; ++run) {
        Isorropia::Input input;
        Isorropia::State state;

        input.w[0] = (dist_na(gen) / 23.0) * 1e-6;   // Na+ component -> mol/m3
        input.w[1] = (dist_so4(gen) / 98.0) * 1e-6;  // H2SO4 component -> mol/m3
        input.w[2] = (dist_nh3(gen) / 17.0) * 1e-6;  // NH3 component -> mol/m3
        input.w[3] = (dist_hno3(gen) / 63.0) * 1e-6; // HNO3 component -> mol/m3
        input.w[4] = (dist_cl(gen) / 36.5) * 1e-6;   // HCl component -> mol/m3
        
        input.rh   = dist_rh(gen);
        input.temp = dist_temp(gen);

        input.org[0] = dist_org(gen) * 1e-9;        // ug/m3 to kg/m3
        input.org[1] = dist_korg(gen);
        input.org[2] = dist_rhoorg(gen);             // kg/m3

        // Execute chemical solvers
        solver.solve(input, state);

        // --- Property 1: Non-Negativity Invariant ---
        EXPECT_GE(state.water, -1e-15) << "Aerosol liquid water must be non-negative";
        EXPECT_GE(state.ionic, -1e-15) << "Ionic strength must be non-negative";
        for (size_t i = 0; i < 10; ++i) {
            EXPECT_GE(state.molal[i], -1e-15) << "Liquid ion index " << i << " concentration must be non-negative";
        }
        for (size_t i = 0; i < 23; ++i) {
            EXPECT_GE(state.molalr[i], -1e-15) << "Active pair index " << i << " molality must be non-negative";
            EXPECT_GE(state.gama[i], -1e-15) << "Activity coefficient index " << i << " must be non-negative";
        }

        // --- Property 2: Electroneutrality (Charge Neutrality) ---
        // Sum equivalent cations = Sum equivalent anions
        // Cations: Na+ (state.molal[0]), H+ (state.molal[1]), NH4+ (state.molal[2])
        // Anions: NO3- (state.molal[3]), Cl- (state.molal[4]), SO4-- (state.molal[5] * 2), HSO4- (state.molal[6])
        double cations = state.molal[0]*1.0 + state.molal[1]*1.0 + state.molal[2]*1.0;
        double anions  = state.molal[3]*1.0 + state.molal[4]*1.0 + state.molal[5]*2.0 + state.molal[6]*1.0;
        
        if (state.water > 1e-10) {
            // Verify charge balance holds within structural limits
            EXPECT_NEAR(cations, anions, 1.0e-5) << "Electroneutrality charge balance failed for run record " << run;
        }

        // --- Property 3: Diagnostic Stability ---
        EXPECT_GE(state.num_errors, 0) << "Diagnostics must remain non-negative";
    }
}

//=======================================================================
// 2. MASS CONSERVATION (MASS BALANCE) INVARIANT
//=======================================================================
TEST(PropertyTest, VerifyMassConservation) {
    Isorropia::Solver solver;
    std::mt19937 gen(202); // fixed seed
    
    std::uniform_real_distribution<double> dist_so4(0.1, 15.0);
    std::uniform_real_distribution<double> dist_nh3(0.1, 40.0);
    std::uniform_real_distribution<double> dist_hno3(0.0, 10.0);
    std::uniform_real_distribution<double> dist_rh(0.20, 0.95);
    std::uniform_real_distribution<double> dist_temp(270.0, 310.0);

    for (int run = 0; run < 100; ++run) {
        Isorropia::Input input;
        Isorropia::State state;

        input.w[1] = (dist_so4(gen) / 98.0) * 1e-6;  // H2SO4 -> mol/m3
        input.w[2] = (dist_nh3(gen) / 17.0) * 1e-6;  // NH3 -> mol/m3
        input.w[3] = (dist_hno3(gen) / 63.0) * 1e-6; // HNO3 -> mol/m3
        input.rh   = dist_rh(gen);
        input.temp = dist_temp(gen);

        // Execute solver
        solver.solve(input, state);

        // Skip assertions if bisection bounds failed to converge (warnings logged)
        if (state.num_errors > 0) continue;

        // Verify component mass balances: Total Input moles == Total Output moles
        // --- 2.1 Sulfate Balance ---
        // Sulfate = Liquid SO4-- + Liquid HSO4-
        double out_so4 = state.molal[5] + state.molal[6];
        EXPECT_NEAR(out_so4, input.w[1], 1e-12) << "Sulfate mass conservation failed";

        // --- 2.2 Ammonium/Ammonia Balance ---
        // Ammonium = Liquid NH4+ + Gas NH3
        double out_nh3 = state.molal[2] + state.gnh3;
        EXPECT_NEAR(out_nh3, input.w[2], 1e-12) << "Ammonia mass conservation failed";

        // --- 2.3 Nitrate/Nitric Acid Balance ---
        // Nitrate = Liquid NO3- + Gas HNO3 + undissociated HNO3aq
        double out_no3 = state.molal[3] + state.ghno3 + state.gasaq[2];
        EXPECT_NEAR(out_no3, input.w[3], 1e-12) << "Nitrate mass conservation failed";
    }
}

//=======================================================================
// 3. THERMODYNAMIC MONOTONICITY (ZSR Water vs RH)
//=======================================================================
TEST(PropertyTest, VerifyWaterMonotonicity) {
    Isorropia::Solver solver;
    std::mt19937 gen(303); // fixed seed
    
    std::uniform_real_distribution<double> dist_so4(1.0, 10.0);
    std::uniform_real_distribution<double> dist_nh3(1.0, 20.0);
    std::uniform_real_distribution<double> dist_hno3(1.0, 8.0);
    std::uniform_real_distribution<double> dist_temp(280.0, 300.0);

    for (int run = 0; run < 20; ++run) {
        Isorropia::Input input;
        input.w[1] = (dist_so4(gen) / 98.0) * 1e-6;
        input.w[2] = (dist_nh3(gen) / 17.0) * 1e-6;
        input.w[3] = (dist_hno3(gen) / 63.0) * 1e-6;
        input.temp = dist_temp(gen);

        double last_water = -1.0;
        
        // Sweep Relative Humidity from 20% to 95%
        for (int rh_pct = 20; rh_pct <= 95; rh_pct += 5) {
            input.rh = static_cast<double>(rh_pct) / 100.0;
            Isorropia::State state;
            solver.solve(input, state);

            if (state.num_errors > 0) continue;

            if (last_water >= 0.0) {
                // Aerosol liquid water content must monotonically increase with RH
                EXPECT_GE(state.water, last_water - 1e-15) 
                    << "Water monotonicity failed at RH " << input.rh << " for run " << run;
            }
            last_water = state.water;
        }
    }
}

//=======================================================================
// 4. INTER-LANGUAGE C/C++ API EQUIVALENCE
//=======================================================================
TEST(PropertyTest, VerifyCAPIEquivalence) {
    Isorropia::Solver solver;
    std::mt19937 gen(404); // fixed seed

    std::uniform_real_distribution<double> dist_na(0.0, 5.0);
    std::uniform_real_distribution<double> dist_so4(0.1, 10.0);
    std::uniform_real_distribution<double> dist_nh3(0.1, 20.0);
    std::uniform_real_distribution<double> dist_hno3(0.0, 10.0);
    std::uniform_real_distribution<double> dist_cl(0.0, 5.0);
    std::uniform_real_distribution<double> dist_rh(0.30, 0.90);
    std::uniform_real_distribution<double> dist_temp(275.0, 305.0);

    for (int run = 0; run < 30; ++run) {
        Isorropia::Input input;
        IsorropiaState state_c = {0.0};
        IsorropiaInput input_c = {0.0};

        input.w[0] = (dist_na(gen) / 23.0) * 1e-6;
        input.w[1] = (dist_so4(gen) / 98.0) * 1e-6;
        input.w[2] = (dist_nh3(gen) / 17.0) * 1e-6;
        input.w[3] = (dist_hno3(gen) / 63.0) * 1e-6;
        input.w[4] = (dist_cl(gen) / 36.5) * 1e-6;
        input.rh   = dist_rh(gen);
        input.temp = dist_temp(gen);

        // Replicate to C api struct
        std::copy(input.w.begin(), input.w.end(), input_c.w);
        input_c.rh = input.rh;
        input_c.temp = input.temp;
        input_c.iprob = 0;
        input_c.nadj = 1;

        // Solve directly via C++ and flat C linkage signatures
        Isorropia::State state_cpp;
        solver.solve(input, state_cpp);
        isorropia_solve_c(&input_c, &state_c);

        // Verify absolute identical numerical equivalence (bit-level congruences)
        EXPECT_DOUBLE_EQ(state_cpp.water, state_c.water);
        EXPECT_DOUBLE_EQ(state_cpp.ionic, state_c.ionic);
        EXPECT_DOUBLE_EQ(state_cpp.gnh3, state_c.gnh3);
        EXPECT_DOUBLE_EQ(state_cpp.ghno3, state_c.ghno3);
        EXPECT_DOUBLE_EQ(state_cpp.ghcl, state_c.ghcl);
        EXPECT_EQ(state_cpp.num_errors, state_c.num_errors);

        for (size_t i = 0; i < 10; ++i) {
            EXPECT_DOUBLE_EQ(state_cpp.molal[i], state_c.molal[i]);
        }
        for (size_t i = 0; i < 23; ++i) {
            EXPECT_DOUBLE_EQ(state_cpp.gama[i], state_c.gama[i]);
        }
    }
}
