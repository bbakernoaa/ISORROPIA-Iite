# ISORROPIA-Lite C++ Port (Phase 3: Activity Coefficients, ZSR Iterations, and Speciation Solvers) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the dynamic Kusik-Meissner activity coefficient models (`KMTAB` and tables), Bromley equations, iterative multi-component ZSR water calculators (`CALCMR`), and the mathematical speciation loops for Case 1 metastables to transition from mocked calculations to full physics.

**Architecture:** Compile activities, Bromley formulations, and ZSR models in standard modules. All calculations are encapsulated inside the local `Isorropia::State` to prevent global memory pollution.

**Tech Stack:** C++17, CMake, GoogleTest, Python 3

## Global Constraints

- Numerical Equivalence: Tolerance-based equality ($10^{-12}$ to $10^{-15}$).
- Thread Safety: No global variables or COMMON blocks.
- Robustness: Zero usage of C++ exceptions (`throw`) in solvers.
- Library Structure: Standalone C++ Compiled Library (C++17) using CMake.
- 0-Based Indexing for all C++ arrays.
- Scientist Readability: All variables, units, and chemical equations must be extensively documented in comments mapping directly back to original Fortran names and the scientific literature.

---

### Task 1: Tabulated Kusik-Meissner Activity Coefficients

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare activity structures)
- Create: `src/ActivityCoefficients.cpp` (implement KMTAB and pre-tabulated grids)
- Create: `tests/activity_parser.py` (optional automation script for grid loading)
- Modify: `CMakeLists.txt` (compile ActivityCoefficients.cpp)
- Create: `tests/test_activities.cpp` (test coefficient selections)
- Modify: `tests/CMakeLists.txt` (register tests)

**Interfaces:**
- Consumes: `Isorropia::State::ionic` (ionic strength) and `Isorropia::State::temp` (temperature).
- Produces: `void State::km_tab(double ionic_strength, double temp_k, std::array<double, 23>& g0)` setting binary activity coefficients.

- [ ] **Step 1: Declare km_tab inside Solver.hpp**

```cpp
    /**
     * @brief Computes binary activity coefficients using pre-tabulated Kusik-Meissner grids.
     * 
     * Replaces Fortran 'SUBROUTINE KMTAB'.
     */
    void km_tab(double ionic_strength, double temp_k, std::array<double, 23>& g0);
```

- [ ] **Step 2: Write tests/activity_parser.py to extract data from Fortran**

Replicate the pre-tabulated coefficient grids for `KM198`, `KM223`, `KM248`, `KM273`, `KM298`, and `KM323` inside `src/ActivityCoefficients.cpp` using an extraction script or direct static mappings.

- [ ] **Step 3: Implement km_tab in src/ActivityCoefficients.cpp**

```cpp
#include "Isorropia/Solver.hpp"
#include <cmath>
#include <algorithm>

namespace Isorropia {

void State::km_tab(double ionic_strength, double temp_k, std::array<double, 23>& g0) {
    // Determine temperature range index IND
    int ind = static_cast<int>((temp_k - 198.0) / 25.0 + 0.5) + 1;
    ind = std::max(1, std::min(ind, 6));

    // Call appropriate pre-tabulated routine and populate g0 array
    // ...
}

} // namespace Isorropia
```

- [ ] **Step 4: Write tests in tests/test_activities.cpp**

```cpp
#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(ActivityTest, KusikMeissnerAt298) {
    Isorropia::State state;
    std::array<double, 23> g0 = {0.0};
    state.km_tab(4.295, 298.15, g0);
    EXPECT_GT(g0[3], 0.0); // Check non-zero coefficients
}
```

- [ ] **Step 5: Compile and verify unit tests**

Run: `cd build && cmake .. && make && ctest -V`
Expected: Activity assertions PASS.

- [ ] **Step 6: Commit**

```bash
git add include/Isorropia/Solver.hpp src/ActivityCoefficients.cpp tests/test_activities.cpp CMakeLists.txt tests/CMakeLists.txt
git commit -m "feat: implement tabulated Kusik-Meissner activity coefficient models"
```

---

### Task 2: Multicomponent Activity Solver (CALCACT1)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare activity methods)
- Modify: `src/ActivityCoefficients.cpp` (implement cal_act1 multicomponent calculations)
- Modify: `tests/test_activities.cpp` (add unit tests verifying multi-species activities)

**Interfaces:**
- Consumes: `State::molal`, `State::water`, `State::temp`
- Produces: `void State::cal_act1()` populating multicomponent activity array `State::gama`.

- [ ] **Step 1: Declare cal_act1 in Solver.hpp**

```cpp
    /**
     * @brief Computes multicomponent activity coefficients for Case 1.
     * 
     * Replaces Fortran 'SUBROUTINE CALCACT1'.
     */
    void cal_act1();
```

- [ ] **Step 2: Implement cal_act1 inside src/ActivityCoefficients.cpp**

Port multicomponent Bromley / Kusik-Meissner mix formulations (Fountoukis & Nenes, 2007) from `isocom.f` line 4517:
- Calculate ionic strength of solution.
- Call `km_tab` to fetch binary activity coefficients.
- Apply Debye-Huckel constants $A_\gamma = 0.511 \times (298.15/T)^{1.5}$.
- Update multicomponent activity array `gama`.

- [ ] **Step 3: Add unit tests to tests/test_activities.cpp**

```cpp
TEST(ActivityTest, MulticomponentActivity1) {
    Isorropia::State state;
    state.water = 8.088;
    state.temp = 298.15;
    state.molal[1] = 2.170e-5; // H+
    state.molal[2] = 2.381e-2; // NH4+
    state.molal[5] = 1.014e-2; // SO4--
    
    state.cal_act1();
    
    EXPECT_GT(state.gama[6], 0.0); // gama of active pair
}
```

- [ ] **Step 4: Verify build and test**

Run: `cd build && make && ctest -V`
Expected: Passes multicomponent activity validation.

- [ ] **Step 5: Commit**

```bash
git add include/Isorropia/Solver.hpp src/ActivityCoefficients.cpp tests/test_activities.cpp
git commit -m "feat: implement multicomponent activity coefficient solver"
```

---

### Task 3: Iterative Water Activity & ZSR Iterations (CALCMR)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare water activity calculations)
- Modify: `src/Thermodynamics.cpp` (implement cal_cmr iterations)
- Modify: `tests/test_thermo.cpp` (add ZSR unit tests)

**Interfaces:**
- Consumes: `State::molalr` (molarity of active pairs), meteorological copies.
- Produces: `void State::cal_cmr()` updating total liquid water content `State::water`.

- [ ] **Step 1: Declare cal_cmr inside Solver.hpp**

```cpp
    /**
     * @brief Iteratively solves multi-component ZSR water activity relations.
     * 
     * Replaces Fortran 'SUBROUTINE CALCMR'.
     */
    void cal_cmr();
```

- [ ] **Step 2: Implement cal_cmr inside src/Thermodynamics.cpp**

Port lines 3565-3850 from `isocom.f` replicating ZSR liquid aerosol water calculations. Include organic water contributions `watcmp[23]` using hygroscopicity relations.

- [ ] **Step 3: Write ZSR unit tests in tests/test_thermo.cpp**

```cpp
TEST(ThermodynamicsTest, ZSRWaterUptake) {
    Isorropia::State state;
    state.rh = 0.80;
    state.temp = 298.15;
    state.molalr[3] = 1.750; // (NH4)2SO4 concentration
    
    state.cal_cmr();
    EXPECT_GT(state.water, 0.0); // water uptake must be positive
}
```

- [ ] **Step 4: Run tests**

Run: `cd build && make && ctest -V`
Expected: Water uptake assertions PASS.

- [ ] **Step 5: Commit**

```bash
git add include/Isorropia/Solver.hpp src/Thermodynamics.cpp tests/test_thermo.cpp
git commit -m "feat: implement iterative ZSR water activity and organic water model"
```

---

### Task 4: Metasable Pure Case 1 Solver SPECIATION

**Files:**
- Modify: `src/ForwardSolvers.cpp` (implement analytical and bisection loops of isrp1f)
- Modify: `tests/test_solver.cpp` (unit tests verifying raw mathematical solvers)

**Interfaces:**
- Consumes: `Isorropia::Input`
- Produces: Numerically computed concentrations for NH4-SO4 systems.

- [ ] **Step 1: Write pure math speciation solvers inside src/ForwardSolvers.cpp**

Replace the skeleton/mock values inside `isrp1f` and `isrp2f` with the complete numerical solvers (e.g., implementing bisulfate dissociation, analytical hydrogen activity loops, and bisection roots).

- [ ] **Step 2: Write tests checking mathematical convergence in tests/test_solver.cpp**

Test with varying RH (e.g. 50% to 90% RH) and assert that solutions converge without errors.

- [ ] **Step 3: Verify the whole test suite and run E2E python regression**

Run: `cd build && cmake .. && make && ctest -V && python3 ../tests/regression_runner.py`
Expected: All tests and the regression validation harness PASS successfully.

- [ ] **Step 4: Commit**

```bash
git add src/ForwardSolvers.cpp tests/test_solver.cpp
git commit -m "feat: replace mock solvers with high-fidelity speciation solvers"
```
