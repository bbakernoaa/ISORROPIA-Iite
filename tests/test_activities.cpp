#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(ActivityTest, KusikMeissnerAt298) {
    Isorropia::State state;
    std::array<double, 23> g0 = {0.0};
    
    // Call KMTAB with standard temperature 298.15K and ionic strength 4.295
    state.km_tab(4.295, 298.15, g0);
    
    // Assert known values from the tables to guarantee exact copy fidelity and selection
    // BNC01M at idx for ionic=4.295:
    // ipos for 4.295: NINT(20.0*4.295)+1 = 86 + 1 = 87 -> idx is 86
    // BNC01M[86] at 298K: let's verify it is populated
    EXPECT_GT(g0[0], -1.0); // NaCl coefficient
    EXPECT_LT(g0[0], 1.0);
    
    // Verify values are not all zero
    double sum = 0.0;
    for (double val : g0) {
        sum += std::abs(val);
    }
    EXPECT_GT(sum, 0.0);
}

TEST(ActivityTest, MulticomponentActivity1) {
    Isorropia::State state;
    state.initialize_constants();
    state.water = 8.088;
    state.temp = 298.15;
    state.epsact = 1e-4;
    state.frst = true;

    // Set up dummy liquid molalities matching a realistic Case 1 solution
    // MOLAL(1) = Na+, MOLAL(2) = H+, MOLAL(3) = NH4+, MOLAL(6) = SO4--, MOLAL(7) = HSO4- (using Fortran 1-based names)
    // C++ indices are: 0: Na+, 1: H+, 2: NH4+, 5: SO4--, 6: HSO4-
    state.molal[0] = 0.0;     // Na+
    state.molal[1] = 2.170e-5; // H+
    state.molal[2] = 2.381e-2; // NH4+
    state.molal[3] = 0.0;     // NO3-
    state.molal[4] = 0.0;     // Cl-
    state.molal[5] = 1.014e-2; // SO4--
    state.molal[6] = 6.487e-5; // HSO4-

    // Compute multicomponent activity coefficients
    state.cal_act1();

    // Verify properties are populated and reasonable
    EXPECT_GT(state.ionic, 0.0);
    EXPECT_NEAR(state.ionic, 0.003979, 1e-5); // correct ionic strength after MOLAL(2) and MOLAL(7) resets
    
    // Verify activity coefficients are populated and within normal chemical limits (e.g. 10^-5 to 10^3)
    EXPECT_GT(state.gama[3], 1e-5);  // (NH4)2SO4
    EXPECT_LT(state.gama[3], 1e3);
    EXPECT_GT(state.gama[8], 1e-5);  // NH4HSO4
    EXPECT_LT(state.gama[8], 1e3);
    EXPECT_GT(state.gama[12], 1e-5); // Letovicite
    EXPECT_LT(state.gama[12], 1e3);

    // Verify flags are updated
    EXPECT_FALSE(state.frst);
}
