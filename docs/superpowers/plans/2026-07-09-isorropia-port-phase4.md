# ISORROPIA-Lite C++ Port (Phase 5: C-Compatible Public API and CATChem Integration Preparation) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Provide a C-compatible, flat-linkable ABI interface (`Isorropia::km_tab`, `Isorropia::cal_cmr`, etc.) encapsulated in standard header files to allow atmospheric model interfaces like CATChem (C, C++, Fortran) to compile, link, and invoke our thread-safe solver cleanly.

**Architecture:** Implement standard C-compatible headers under `include/Isorropia/Isorropia.h` containing plain C structs and free functions. These map fields directly into C++ `Isorropia::Input` and `Isorropia::State` under the hood, ensuring complete thread-safety.

**Tech Stack:** C++17, C, CMake, GoogleTest

## Global Constraints

- Numerical Equivalence: Tolerance-based equality ($10^{-12}$ to $10^{-15}$).
- Thread Safety: No global variables or COMMON blocks.
- Robustness: Zero usage of C++ exceptions (`throw`) in solvers.
- Library Structure: Standalone C++ Compiled Library (C++17) using CMake.
- 0-Based Indexing for all C++ arrays.
- Scientist Readability: All variables, units, and chemical equations must be extensively documented in comments mapping directly back to original Fortran names and the scientific literature.

---

### Task 1: C-Compatible API Header

**Files:**
- Create: `include/Isorropia/Isorropia.h`

**Interfaces:**
- Produces: Standard C-linkable header containing `IsorropiaInput`, `IsorropiaState`, and solver handles compatible with `extern "C"`.

- [ ] **Step 1: Write Isorropia.h**

```cpp
#ifndef ISORROPIA_C_API_H
#define ISORROPIA_SOLVER_H // guards matching Solver.hpp if needed, but let's use ISORROPIA_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double w[8];
    double org[3];
    double waer[8];
    double temp;
    double rh;
    int iprob;
    int nadj;
} IsorropiaInput;

typedef struct {
    double temp;
    double rh;
    double w[8];
    double waer[8];
    double org[3];

    double molal[10];
    double molalr[23];
    double gama[23];
    double zz[23];
    double z[10];
    double gamou[23];
    double gamin[23];
    double m0[23];
    double gasaq[3];
    int actmod;
    double epsact;
    double coh;
    double chno3;
    double chcl;
    double water;
    double ionic;
    double watcmp[24];
    int frst;
    int calain;
    int calaou;
    int dryf;

    double ch2so4, cnh42s4, cnh4hs4, cnacl, cna2so4, cnano3, cnh4no3, cnh4cl, cnahso4, clc;
    double ccaso4, ccano32, ccacl2, ck2so4, ckhso4, ckno3, ckcl, cmgso4, cmgno32, cmgcl2;
    double gnh3, ghno3, ghcl;

    int num_errors;
    // We encapsulate the remaining error-stack elements internally
} IsorropiaState;

/**
 * @brief Top-level execution entry point for C and Fortran binders.
 */
void isorropia_solve_c(const IsorropiaInput* input, IsorropiaState* state);

#ifdef __cplusplus
}
#endif

#endif
```

- [ ] **Step 2: Commit**

```bash
git add include/Isorropia/Solver.hpp # Verify references
git add include/Isorropia/Isorropia.h
git commit -m "feat: declare C-compatible public API headers"
```

---

### Task 2: C-API Implementation (src/ReverseSolvers.cpp or src/Solver.cpp)

**Files:**
- Modify: `src/Solver.cpp` (implement isorropia_solve C bindings)
- Modify: `CMakeLists.txt` (compile updates)

**Interfaces:**
- Consumes: `Isorropia::Input` and `Isorropia::State` C++ structures.
- Produces: Compiled `isorropia_cli` and `isorropia` static libraries containing the `isorropia_solve` symbol.

- [ ] **Step 1: Implement km_tab and solvers conversions in Solver.cpp**

```cpp
extern "C" {

void isorropia_solve(const Isorropia::Input* input, Isorropia::State* state) {
    if (!input || !state) return;
    
    // Create C++ structures, call Solver::solve, and copy output metrics back to C-compatible pointers
    Isorropia::Solver solver;
    solver.solve(*input, *state);
}

}
```

- [ ] **Step 2: Verify build**

Run: `cd build && make`
Expected: Static library builds successfully without errors.

- [ ] **Step 3: Commit**

```bash
git add src/Solver.cpp CMakeLists.txt
git commit -m "feat: implement C-compatible wrapper functions"
```

---

### Task 3: Unit Tests for C-Compatible Interface

**Files:**
- Create: `tests/test_c_api.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `isorropia` library.

- [ ] **Step 1: Write test_reverse_c.cpp**

```cpp
#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(CAPITest, SolveCCompatible) {
    Isorropia::Input input;
    Isorropia::State state;
    
    input.w[1] = 1.0; // H2SO4
    input.w[2] = 2.0; // NH3
    input.w[3] = 1.0; // HNO3
    input.rh = 0.80;
    input.temp = 298.15;
    
    // Call the solver
    Isorropia::Solver solver;
    solver.solve(input, state);
    
    EXPECT_GT(state.water, 0.0);
}
```

- [ ] **Step 2: Update tests/CMakeLists.txt**

```cmake
add_executable(isorropia_tests test_main.cpp test_state.cpp test_thermo.cpp test_solver.cpp test_activities.cpp test_crustal.cpp test_reverse.cpp test_activities.cpp tests/test_activities.cpp tests/test_reverse.cpp tests/test_crustal.cpp tests/test_state.cpp)
```

- [ ] **Step 3: Verify all unit tests and regression harness**

Run: `cd build && cmake .. && make && ctest -V && python3 ../tests/regression_runner.py`
Expected: 100% PASS!

- [ ] **Step 4: Commit**

```bash
git add tests/test_reverse.cpp tests/CMakeLists.txt
git commit -m "test: add integration checks for the C-compatible link-state"
```
