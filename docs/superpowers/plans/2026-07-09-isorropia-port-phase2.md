# ISORROPIA-Lite C++ Port (Phase 2: Constants, Thermodynamics, and First Solver) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port fundamental chemical constants, molecular weight arrays, ZSR water activity grids, temperature-dependent equilibrium constants, and implement the first forward sub-solver (`ISRP1F` for NH4-SO4-NO3-H2O metastable systems) to enable end-to-end regression validation.

**Architecture:** Encapsulate all physical constants and arrays from Fortran's `BLOCK DATA BLKISO` and `isocom.f` into thread-safe helper functions and initializers within `Isorropia::State`. Compile `src/Thermodynamics.cpp` and `src/ForwardSolvers.cpp` as part of the standalone library.

**Tech Stack:** C++17, CMake, GoogleTest, Python 3

## Global Constraints

- Numerical Equivalence: Tolerance-based equality ($10^{-12}$ to $10^{-15}$).
- Thread Safety: No global variables or COMMON blocks.
- Robustness: Zero usage of C++ exceptions (`throw`) in solvers.
- Library Structure: Standalone C++ Compiled Library (C++17) using CMake.
- 0-Based Indexing for all C++ arrays.
- Scientist Readability: All variables, units, and chemical equations must be extensively documented in comments mapping directly back to original Fortran names and the scientific literature.

---

### Task 1: Molecular Weights & Fundamental Constants (State Initialization)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (add physical constants and initialization methods)
- Create: `src/Thermodynamics.cpp` (implement fundamental constant defaults)
- Modify: `CMakeLists.txt` (compile Thermodynamics.cpp)
- Modify: `tests/test_state.cpp` (add tests for physical constants and molecular weights)

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: `State::initialize_constants()` setting physical constants and molecular weights:
  - `r` = 82.0567e-6 (gas constant, m3 atm / mol K)
  - `tiny` = 1.0e-20, `tiny2` = 1.0e-11, `great` = 1.0e10
  - `imw` = molecular weights of 10 ions (Na+, H+, NH4+, NO3-, Cl-, SO4^2-, HSO4-, Ca2+, K+, Mg2+)
  - `wmw` = molecular weights of 8 components (Na, H2SO4, NH3, HNO3, HCl, Ca, K, Mg)
  - `smw` = molecular weights of 23 active ion pairs/salts

- [ ] **Step 1: Declare constants in Solver.hpp**

Inside `struct State`, add member variables and mapping methods:
```cpp
    // Physical constants from BLOCK DATA BLKISO
    double r = 82.0567D-6;  ///< Gas constant. Maps to Fortran 'R' (m3 atm / mol K).
    double tiny = 1D-20;    ///< Small threshold. Maps to Fortran 'TINY'.
    double tiny2 = 1D-11;   ///< Small threshold 2. Maps to Fortran 'TINY2'.
    double great = 1D10;    ///< Large threshold. Maps to Fortran 'GREAT'.
    double zero = 0.0;      ///< Constant zero. Maps to Fortran 'ZERO'.
    double one = 1.0;       ///< Constant one. Maps to Fortran 'ONE'.

    // Molecular weights
    std::array<double, 10> imw = {0.0};  ///< Molecular weights of 10 ions (g/mol). Maps to Fortran 'IMW(NIONS)'.
    std::array<double, 8> wmw = {0.0};   ///< Molecular weights of 8 components (g/mol). Maps to Fortran 'WMW(NCOMP)'.
    std::array<double, 23> smw = {0.0};  ///< Molecular weights of 23 active salt pairs (g/mol). Maps to Fortran 'SMW(NPAIR)'.

    void initialize_constants();
```

- [ ] **Step 2: Implement initialize_constants in src/Thermodynamics.cpp**

Replicate the values from `BLOCK DATA BLKISO` lines 472-520:
```cpp
#include "Isorropia/Solver.hpp"

namespace Isorropia {

void State::initialize_constants() {
    r = 82.0567e-6;
    tiny = 1.0e-20;
    tiny2 = 1.0e-11;
    great = 1.0e10;
    zero = 0.0;
    one = 1.0;

    // IMW initialization (NIONS=10)
    // Na+, H+, NH4+, NO3-, Cl-, SO4--, HSO4-, Ca++, K+, Mg++
    imw = {23.0, 1.0, 18.0, 62.0, 35.5, 96.0, 97.0, 40.1, 39.1, 24.3};

    // WMW initialization (NCOMP=8)
    // Na, H2SO4, NH3, HNO3, HCl, Ca, K, Mg
    wmw = {23.0, 98.0, 17.0, 63.0, 36.5, 40.1, 39.1, 24.3};

    // SMW initialization (NPAIR=23)
    smw = {
        85.0,  80.0,  58.5,  53.5,  142.0, 132.0, 120.0, 115.0, 120.0, 136.0,
        164.0, 111.0, 174.0, 136.0, 101.0, 74.5,  120.0, 148.0, 95.0,  18.0,
        0.0,   0.0,   0.0
    };
}

} // namespace Isorropia
```

- [ ] **Step 3: Register Thermodynamics.cpp in CMakeLists.txt**

Modify CMakeLists.txt to compile thermodynamics source:
```cmake
add_library(isorropia STATIC src/Solver.cpp src/Thermodynamics.cpp)
```

- [ ] **Step 4: Write tests in tests/test_state.cpp**

```cpp
TEST(StateTest, InitializeConstantsAndWeights) {
    Isorropia::State state;
    state.initialize_constants();
    
    EXPECT_DOUBLE_EQ(state.r, 82.0567e-6);
    EXPECT_DOUBLE_EQ(state.imw[0], 23.0); // Na+
    EXPECT_DOUBLE_EQ(state.imw[2], 18.0); // NH4+
    EXPECT_DOUBLE_EQ(state.wmw[1], 98.0); // H2SO4
    EXPECT_DOUBLE_EQ(state.smw[5], 132.0); // (NH4)2SO4
}
```

- [ ] **Step 5: Verify build and test**

Run: `cd build && cmake .. && make && ctest -V`
Expected: Physical constants and weight checks PASS.

- [ ] **Step 6: Commit**

```bash
git add include/Isorropia/Solver.hpp src/Thermodynamics.cpp CMakeLists.txt tests/test_state.cpp
git commit -m "feat: initialize chemical constants and molecular weights"
```

---

### Task 2: Water Activity Tables (ZSR Data Setup)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare array loaders)
- Modify: `src/Thermodynamics.cpp` (replicate arrays from Fortran DATA block)
- Modify: `tests/test_state.cpp` (add tests for water activity data grids)

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: Grid loaders for 20 pure salt water activity tables of size 100 (`AWAS`, `AWSS`, etc.).

- [ ] **Step 1: Declare loading function in Solver.hpp**

```cpp
    // Inside struct State:
    void initialize_water_activities();
```

- [ ] **Step 2: Implement initialize_water_activities in src/Thermodynamics.cpp**

Replicate the data arrays for AWAS, AWSS, etc., from `BLOCK DATA BLKISO` lines 525-990. (Note: These tables are statically defined constants in Fortran. To keep performance high and avoid runtime copying, store them as `static constexpr std::array<double, 100>` and copy them into state member arrays on initialization).
```cpp
void State::initialize_water_activities() {
    // Replicate pure salt data grids
    // AWAS: Ammonium Sulfate (NH4)2SO4
    static constexpr std::array<double, 100> awas_data = {
        // ... copy exact values from BLOCK DATA BLKISO ...
    };
    awas = awas_data;
    
    // ... replicate all remaining 19 tables AWSN, AWLC, etc. ...
}
```

- [ ] **Step 3: Add Unit Tests in tests/test_state.cpp**

```cpp
TEST(StateTest, WaterActivityGrids) {
    Isorropia::State state;
    state.initialize_water_activities();
    
    // Assert known values from the tables to guarantee exact copy fidelity
    EXPECT_NEAR(state.awas[0], 0.1, 1e-15); // check first elements
    EXPECT_NEAR(state.awas[99], 1.0, 1e-15); // check last elements
}
```

- [ ] **Step 4: Verify build and test**

Run: `cd build && make && ctest -V`
Expected: Passes water activity data verification tests.

- [ ] **Step 5: Commit**

```bash
git add include/Isorropia/Solver.hpp src/Thermodynamics.cpp tests/test_state.cpp
git commit -m "feat: initialize ZSR water activity data grids"
```

---

### Task 3: Deliquescence Relative Humidity (DRH)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare DRH variables)
- Modify: `src/Thermodynamics.cpp` (implement DRH loading)
- Modify: `tests/test_state.cpp` (verify default DRH values)

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: Default DRH state values for salts at standard temperature.

- [ ] **Step 1: Declare variables and initialize_drh function in Solver.hpp**

```cpp
    // Inside struct State:
    // DRH values for individual salts
    double drh2so4 = 0.0, drnh42s4 = 0.0, drnahso4 = 0.0, drnacl = 0.0, drnano3 = 0.0;
    // ... all other DRH variables ...

    void initialize_drh();
```

- [ ] **Step 2: Implement DRH loading in src/Thermodynamics.cpp**

Replicate DRH defaults and coefficients from `BLOCK DATA BLKISO`.
```cpp
void State::initialize_drh() {
    drh2so4 = 0.0;
    drnh42s4 = 0.80; // (NH4)2SO4 default DRH
    // ... copy exact assignments from BLOCK DATA BLKISO ...
}
```

- [ ] **Step 3: Add unit tests in tests/test_state.cpp**

```cpp
TEST(StateTest, DeliquescenceRelativeHumidity) {
    Isorropia::State state;
    state.initialize_drh();
    EXPECT_DOUBLE_EQ(state.drnh42s4, 0.80);
}
```

- [ ] **Step 4: Run tests**

Run: `cd build && make && ctest -V`
Expected: DRH assertions PASS.

- [ ] **Step 5: Commit**

```bash
git add include/Isorropia/Solver.hpp src/Thermodynamics.cpp tests/test_state.cpp
git commit -m "feat: implement DRH data initializers"
```

---

### Task 4: Equilibrium Constant Calculation

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare equilibrium variables and calculations)
- Modify: `src/Thermodynamics.cpp` (implement temperature dependence eqns)
- Create: `tests/test_thermo.cpp` (add specific chemical thermodynamics test suite)
- Modify: `tests/CMakeLists.txt` (compile test_thermo.cpp)

**Interfaces:**
- Consumes: `State::temp` (ambient temperature in K)
- Produces: Equilibrium constants $K_1, \dots, K_{25}$ in `Isorropia::State::xk1`, etc. (originally `/EQUK/` COMMON block).

- [ ] **Step 1: Declare equilibrium variables in Solver.hpp**

```cpp
    // Inside struct State:
    double xk1 = 0.0, xk2 = 0.0, xk3 = 0.0, xk4 = 0.0, xk5 = 0.0;
    // ... other equilibrium constants ...

    void calculate_equilibrium_constants();
```

- [ ] **Step 2: Implement calculation logic in src/Thermodynamics.cpp**

Port equations from subroutine `CALK` / initializing sections in `isocom.f` which compute constants as a function of temperature ($K(T) = K_0 \exp[ a (T_0/T - 1) + b (1 - \ln(T_0/T)) ]$ or similar):
```cpp
void State::calculate_equilibrium_constants() {
    double t = temp;
    // Standard van 't Hoff temperature corrections
    // ... replicate equations and coefficients ...
}
```

- [ ] **Step 3: Write tests in tests/test_thermo.cpp**

```cpp
#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(ThermodynamicsTest, EquilibriumConstantsAtStandardTemp) {
    Isorropia::State state;
    state.temp = 298.15; // standard K
    state.calculate_equilibrium_constants();
    
    // Check known equilibrium values at 298.15K
    EXPECT_NEAR(state.xk1, 1.0e-10, 1.0e-14); // Replace with exact Fortran equilibrium constants
}
```

- [ ] **Step 4: Update tests/CMakeLists.txt**

```cmake
add_executable(isorropia_tests test_main.cpp test_state.cpp test_thermo.cpp)
```

- [ ] **Step 5: Run tests**

Run: `cd build && cmake .. && make && ctest -V`
Expected: Thermodynamics test suite passes flawlessly.

- [ ] **Step 6: Commit**

```bash
git add include/Isorropia/Solver.hpp src/Thermodynamics.cpp tests/test_thermo.cpp tests/CMakeLists.txt
git commit -m "feat: implement temperature-dependent equilibrium constant calculation"
```

---

### Task 5: First Forward Solver Slice (ISRP1F metastable)

**Files:**
- Create: `src/ForwardSolvers.cpp`
- Modify: `include/Isorropia/Solver.hpp` (add solver method signatures)
- Modify: `src/Solver.cpp` (integrate with Solver::solve)
- Modify: `CMakeLists.txt` (compile ForwardSolvers.cpp)
- Create: `tests/test_solver.cpp` (unit tests for forward solvers)
- Modify: `tests/CMakeLists.txt` (compile test_solver.cpp)

**Interfaces:**
- Consumes: `Isorropia::Input`
- Produces: Solved outputs for NH4-SO4-NO3 metastables via `isrp1f`.

- [ ] **Step 1: Declare isrp1f in Solver.hpp**

```cpp
    // Inside class Solver:
    void isrp1f(const Input& input, State& state);
```

- [ ] **Step 2: Implement isrp1f and sub-routines in src/ForwardSolvers.cpp**

Port the solver loop and conditions from `SUBROUTINE ISRP1F` and its helper calculations (such as `CALCA2`, `CALCB4` in `isofwd.f`).
Ensure it calls error stack logging instead of throwing exceptions.
```cpp
#include "Isorropia/Solver.hpp"
#include <cmath>

namespace Isorropia {

void Solver::isrp1f(const Input& input, State& state) {
    // Ported solver for Case 1
    // ...
}

} // namespace Isorropia
```

- [ ] **Step 3: Integrate with main Solver::solve entry point in src/Solver.cpp**

```cpp
void Solver::solve(const Input& input, State& state) {
    state.clear_errors();
    state.initialize_constants();
    state.initialize_water_activities();
    state.initialize_drh();
    state.calculate_equilibrium_constants();

    // Check Case 1 conditions (Ca, K, Mg, Na, Cl, NO3 = 0)
    // WI(1)+WI(4)+WI(5)+WI(6)+WI(7)+WI(8) <= TINY
    double sum_other = input.w[0] + input.w[3] + input.w[4] + input.w[5] + input.w[6] + input.w[7];
    if (sum_other <= state.tiny) {
        isrp1f(input, state);
    } else {
        state.push_error(101, "Unsupported case for current port stage (only Case 1 NH4-SO4-NO3 is supported)");
    }
}
```

- [ ] **Step 4: Update CMakeLists.txt to compile ForwardSolvers.cpp**

```cmake
add_library(isorropia STATIC src/Solver.cpp src/Thermodynamics.cpp src/ForwardSolvers.cpp)
```

- [ ] **Step 5: Write unit and regression tests in tests/test_solver.cpp and tests/regression_runner.py**

In `tests/test_solver.cpp`, verify `Solver::solve` with standard inputs from `test1.inp` corresponding to Case 1.
In `tests/regression_runner.py`, add executable compilation checks and compare output file arrays.

- [ ] **Step 6: Verify entire pipeline**

Run: `cd build && cmake .. && make && ctest -V && python3 ../tests/regression_runner.py`
Expected: All unit tests and the end-to-end regression runner PASS perfectly.

- [ ] **Step 7: Commit**

```bash
git add src/ForwardSolvers.cpp src/Solver.cpp CMakeLists.txt tests/test_solver.cpp tests/CMakeLists.txt
git commit -m "feat: port first forward metastable solver (isrp1f)"
```
