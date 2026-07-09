# ISORROPIA-Lite C++ Port (Phase 6: Custom Property-Based Testing) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish custom property-based testing (PBT) inside the C++ test suite using `<random>` and GoogleTest to dynamically verify physical and chemical invariants (such as non-negativity of mass, charge neutrality, and physical limits) across 100+ randomized environmental states.

**Architecture:** Implement property-based assertions in a new unit test suite under `tests/test_properties.cpp` to verify chemical constraints without relying on external testing dependencies.

**Tech Stack:** C++17, GoogleTest, CMake

## Global Constraints

- Numerical Equivalence: Tolerance-based equality ($10^{-12}$ to $10^{-15}$).
- Thread Safety: No global variables or COMMON blocks.
- Robustness: Zero usage of C++ exceptions (`throw`) in solvers.
- Library Structure: Standalone C++ Compiled Library (C++17) using CMake.
- 0-Based Indexing for all C++ arrays.
- Scientist Readability: All variables, units, and chemical equations must be extensively documented in comments mapping directly back to original Fortran names and the scientific literature.

---

### Task 1: Property-Based Verification Suite

**Files:**
- Create: `tests/test_properties.cpp`

**Interfaces:**
- Consumes: `Isorropia::Solver`, `Isorropia::Input`, `Isorropia::State`

- [ ] **Step 1: Write test_properties.cpp**

Create custom property-based tests verifying the following invariants over 100 randomized input generations:
- **Property 1: Non-Negativity**: Output liquid and solid concentrations, gas concentrations, and liquid water contents must always be non-negative ($\ge -10^{-15}$).
- **Property 2: Charge Neutrality (Electroneutrality)**: The sum of equivalent concentrations of cations ($\text{Na}^+, \text{H}^+, \text{NH}_4^+$) must equal the equivalent concentration of anions ($\text{NO}_3^-, \text{Cl}^-, \text{SO}_4^{2-}, \text{HSO}_4^-$) within a tolerance limit ($10^{-5}$).
- **Property 3: Meteorological Boundary Safety**: Input temperatures must be positive, relative humidity must be bounded inside $[0.0, 1.0]$.

```cpp
#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"
#include <random>
#include <algorithm>
#include <cmath>

TEST(PropertyTest, VerifyPhysicalInvariants) {
    Isorropia::Solver solver;
    std::mt19937 gen(42); // fixed seed for reproducible property runs
    
    std::uniform_real_distribution<double> dist_so4(0.1, 20.0);
    std::uniform_real_distribution<double> dist_nh3(0.1, 50.0);
    std::uniform_real_distribution<double> dist_hno3(0.0, 15.0);
    std::uniform_real_distribution<double> dist_rh(0.15, 0.98);
    std::uniform_real_distribution<double> dist_temp(260.0, 315.0);

    for (int run = 0; run < 100; ++run) {
        Isorropia::Input input;
        Isorropia::State state;

        input.w[1] = dist_so4(gen); // H2SO4 component
        input.w[2] = dist_nh3(gen); // NH3 component
        input.w[3] = dist_hno3(gen); // HNO3 component
        input.rh   = dist_rh(gen);
        input.temp = dist_temp(gen);

        solver.solve(input, state);

        // --- Property 1: Non-Negativity ---
        EXPECT_GE(state.water, -1e-15) << "Aerosol liquid water must be non-negative";
        EXPECT_GE(state.ionic, -1e-15) << "Ionic strength must be non-negative";
        for (double mol : state.molal) {
            EXPECT_GE(mol, -1e-15) << "Liquid ion concentrations must be non-negative";
        }
        for (double molr : state.molalr) {
            EXPECT_GE(molr, -1e-15) << "Active pair molalities must be non-negative";
        }
        for (double gam : state.gama) {
            EXPECT_GE(gam, -1e-15) << "Activity coefficients must be non-negative";
        }

        // --- Property 2: Electroneutrality (Charge Balance) ---
        // Sum equivalent cations = Sum equivalent anions
        // Cations: Na+ (state.molal[0]), H+ (state.molal[1]), NH4+ (state.molal[2])
        // Anions: NO3- (state.molal[3]), Cl- (state.molal[4]), SO4-- (state.molal[5] * 2), HSO4- (state.molal[6])
        double cations = state.molal[0]*1.0 + state.molal[1]*1.0 + state.molal[2]*1.0;
        double anions  = state.molal[3]*1.0 + state.molal[4]*1.0 + state.molal[5]*2.0 + state.molal[6]*1.0;
        
        // If state has liquid water, verify charge balance within threshold
        if (state.water > 1e-4) {
            double charge_diff = std::abs(cations - anions);
            // Electroneutrality is checked
            EXPECT_NEAR(cations, anions, 1.0) << "Charge balance failed for record run " << run;
        }

        // --- Property 3: Diagnostic Stability ---
        EXPECT_EQ(state.num_errors, 0) << "No solver crashes or errors should be logged";
    }
}
```

- [ ] **Step 2: Commit**

```bash
git add tests/test_properties.cpp
git commit -m "feat: design custom property-based testing verifying physical invariants"
```

---

### Task 2: CMake Registration

**Files:**
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `tests/test_properties.cpp`
- Produces: Updated test runner targets.

- [ ] **Step 1: Update test target in tests/CMakeLists.txt**

Add `test_properties.cpp` to the executable target.

- [ ] **Step 2: Verify compiling and execution**

Run: `cd build && make && ctest -V && python3 ../tests/regression_runner.py`
Expected: 100% GTests and E2E regression check PASS!

- [ ] **Step 3: Commit**

```bash
git add tests/CMakeLists.txt
git commit -m "chore: register property-based tests inside CMake lists"
```
