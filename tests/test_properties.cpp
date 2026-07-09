#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"
#include <random>
#include <algorithm>
#include <cmath>

TEST(PropertyTest, VerifyPhysicalInvariants) {
    Isorropia::Solver solver;
    std::mt19937 gen(42); // fixed seed for reproducible property-based testing runs
    
    // Boundary ranges covering standard atmospheric parameters
    std::uniform_real_distribution<double> dist_so4(0.1, 20.0);   // Total sulfate: 0.1 to 20 µg/m3
    std::uniform_real_distribution<double> dist_nh3(0.1, 50.0);   // Total ammonia: 0.1 to 50 µg/m3
    std::uniform_real_distribution<double> dist_hno3(0.0, 15.0);  // Total nitrate: 0.0 to 15 µg/m3
    std::uniform_real_distribution<double> dist_rh(0.15, 0.98);   // RH: 15% to 98%
    std::uniform_real_distribution<double> dist_temp(260.0, 315.0); // Temperature: 260K to 315K

    for (int run = 0; run < 100; ++run) {
        Isorropia::Input input;
        Isorropia::State state;

        input.w[1] = dist_so4(gen);  // H2SO4 component
        input.w[2] = dist_nh3(gen);  // NH3 component
        input.w[3] = dist_hno3(gen); // HNO3 component
        input.rh   = dist_rh(gen);
        input.temp = dist_temp(gen);

        // Execute chemical solvers
        solver.solve(input, state);

        // --- Property 1: Non-Negativity Invariant ---
        // Verify that no physical variable decays into negative regions
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
        
        if (state.water > 1e-4) {
            double charge_diff = std::abs(cations - anions);
            // Verify charge balance holds within structural limits
            EXPECT_NEAR(cations, anions, 1.0) << "Electroneutrality charge balance failed for run record " << run;
        }

        // --- Property 3: Diagnostic Stability ---
        // Checks that no solver crash, division by zero, or NaN-escape values occurred during runtime
        EXPECT_EQ(state.num_errors, 0) << "No solver errors or bisection failures should be logged";
    }
}
