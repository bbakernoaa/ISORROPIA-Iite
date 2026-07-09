# ISORROPIA-Lite C++ Port (Phase 4: Crustal Forward and Reverse Solvers) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the remaining crustal forward solvers (`ISRP3F`, `ISRP4F`) and the entire reverse solver pipeline (`isrp1r` through `isrp4r`) to enable 100% mathematical simulation coverage of all atmospheric aerosol regimes (including forward/reverse problems with Na, Ca, Cl, K, Mg, and organics).

**Architecture:** Implement forward crustal and reverse solver subroutines inside `src/ForwardSolvers.cpp` and `src/ReverseSolvers.cpp`. Register them in the main entry point `Solver::solve()` to dispatch based on `input.iprob` (0 = Forward, 1 = Reverse) and component presence.

**Tech Stack:** C++17, CMake, GoogleTest, Python 3

## Global Constraints

- Numerical Equivalence: Tolerance-based equality ($10^{-12}$ to $10^{-15}$).
- Thread Safety: No global variables or COMMON blocks.
- Robustness: Zero usage of C++ exceptions (`throw`) in solvers.
- Library Structure: Standalone C++ Compiled Library (C++17) using CMake.
- 0-Based Indexing for all C++ arrays.
- Scientist Readability: All variables, units, and chemical equations must be extensively documented in comments mapping directly back to original Fortran names and the scientific literature.

---

### Task 1: Crustal Forward Solvers (ISRP3F and ISRP4F)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare isrp3f, isrp4f)
- Modify: `src/ForwardSolvers.cpp` (implement isrp3f, isrp4f and associated speciation loops)
- Modify: `src/Solver.cpp` (integrate isrp3f/isrp4f dispatch)
- Create: `tests/test_crustal.cpp` (unit tests for crustal regimes)
- Modify: `tests/CMakeLists.txt` (register tests)

**Interfaces:**
- Consumes: `Isorropia::Input` with crustal components (Ca, K, Mg, Na, Cl > 0).
- Produces: `void Solver::isrp3f(const Input&, State&)` and `void Solver::isrp4f(const Input&, State&)` setting correctSpeciation and liquid water contents.

- [ ] **Step 1: Declare crustal solvers in Solver.hpp**

```cpp
    /**
     * @brief Forward solver for Na-NH4-SO4-NO3-Cl-H2O systems (Case 3).
     * 
     * Maps to legacy Fortran 'SUBROUTINE ISRP3F'.
     */
    void isrp3f(const Input& input, State& state);

    /**
     * @brief Forward solver for Na-NH4-SO4-NO3-Cl-Ca-K-Mg-H2O crustal systems (Case 4).
     * 
     * Maps to legacy Fortran 'SUBROUTINE ISRP4F'.
     */
    void isrp4f(const Input& input, State& state);
```

- [ ] **Step 2: Implement isrp3f and isrp4f in src/ForwardSolvers.cpp**

Port speciation systems for crustals from `isofwd.f` line 800-2719.
Maintain standard Newton-Raphson/bisection loops, ensuring error stacks capture non-convergence gracefully.

- [ ] **Step 3: Update Solver.cpp dispatch**

```cpp
    // Inside Solver::solve:
    // If Ca, K, or Mg are present -> call isrp4f
    // Else if Na or Cl are present -> call isrp3f
    // Else if Nitrate is present -> call isrp2f
    // Else -> call isrp1f
```

- [ ] **Step 4: Write unit tests in tests/test_crustal.cpp**

```cpp
#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(CrustalTest, ForwardCrustalSolve) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    Isorropia::State state;
    
    // Set up a standard crustal case
    input.w[0] = 1.0; // Na
    input.w[1] = 1.0; // SO4
    input.w[5] = 0.5; // Ca
    input.rh = 0.80;
    
    solver.solve(input, state);
    EXPECT_GT(state.water, 0.0);
}
```

- [ ] **Step 5: Verify compilation and build**

Run: `cd build && cmake .. && make && ctest -V`
Expected: Crustal assertions PASS.

- [ ] **Step 6: Commit**

```bash
git add include/Isorropia/Solver.hpp src/ForwardSolvers.cpp src/Solver.cpp tests/test_crustal.cpp tests/CMakeLists.txt
git commit -m "feat: implement forward crustal solvers ISRP3F and ISRP4F"
```

---

### Task 2: Reverse Solvers Pipeline & Case 1-2 Reverse Solvers (isrp1r, isrp2r)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare reverse solver methods)
- Create: `src/ReverseSolvers.cpp` (implement isrp1r and isrp2r)
- Modify: `src/Solver.cpp` (integrate reverse problem routing)
- Modify: `CMakeLists.txt` (compile ReverseSolvers.cpp)
- Create: `tests/test_reverse.cpp` (unit tests for reverse solvers)
- Modify: `tests/CMakeLists.txt` (register tests)

**Interfaces:**
- Consumes: `Isorropia::Input` (with iprob = 1)
- Produces: `void Solver::isrp1r(const Input&, State&)` and `void Solver::isrp2r(const Input&, State&)` returning aerosol thermodynamic properties.

- [ ] **Step 1: Declare reverse solvers in Solver.hpp**

```cpp
private:
    void isrp1r(const Input& input, State& state);
    void isrp2r(const Input& input, State& state);
```

- [ ] **Step 2: Implement isrp1r and isrp2r inside src/ReverseSolvers.cpp**

Port reverse calculations from `isorev.f` lines 1-1000. Re-map 1-based arrays to 0-based arrays safely.

- [ ] **Step 3: Update Solver.cpp to dispatch reverse problems**

```cpp
    // Inside Solver::solve:
    if (input.iprob == 1) {
        // Dispatch to reverse solvers: isrp1r, isrp2r, isrp3r, isrp4r
    } else {
        // Dispatch to forward solvers
    }
```

- [ ] **Step 4: Write reverse unit tests in tests/test_reverse.cpp**

Verify that inputs with `iprob = 1` are successfully solved by reverse routines.

- [ ] **Step 5: Verify build**

Run: `cd build && make && ctest -V`
Expected: Reverse unit tests PASS.

- [ ] **Step 6: Commit**

```bash
git add include/Isorropia/Solver.hpp src/ReverseSolvers.cpp src/Solver.cpp CMakeLists.txt tests/test_reverse.cpp tests/CMakeLists.txt
git commit -m "feat: implement reverse solver pipeline and solvers isrp1r/isrp2r"
```

---

### Task 3: Case 3 and 4 Reverse Solvers (isrp3r, isrp4r)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare isrp3r, isrp4r)
- Modify: `src/ReverseSolvers.cpp` (implement isrp3r, isrp4r)
- Modify: `tests/test_reverse.cpp` (add crustal reverse tests)

**Interfaces:**
- Consumes: Crustal aerosol concentrations in input (iprob = 1)
- Produces: `void Solver::isrp3r(const Input&, State&)` and `void Solver::isrp4r(const Input&, State&)`

- [ ] **Step 1: Declare solvers inside Solver.hpp**

```cpp
    void isrp3r(const Input& input, State& state);
    void isrp4r(const Input& input, State& state);
```

- [ ] **Step 2: Implement isrp3r and isrp4r inside src/ReverseSolvers.cpp**

Port crustal reverse speciation math from `isorev.f` lines 1000-1761.

- [ ] **Step 3: Add crustal reverse test assertions in tests/test_reverse.cpp**

Validate correct reverse speciation outputs.

- [ ] **Step 4: Run test suite**

Run: `cd build && make && ctest -V`
Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add include/Isorropia/Solver.hpp src/ReverseSolvers.cpp tests/test_reverse.cpp
git commit -m "feat: implement reverse crustal solvers isrp3r and isrp4r"
```

---

### Task 4: Complete E2E Integration and Multi-file Regression Tests

**Files:**
- Modify: `tests/regression_runner.py` (add support for multiple INP files)

**Interfaces:**
- Consumes: `Partitioning_with_organics.INP`, `Reverse_with_organics.INP`, and corresponding reference outputs.

- [ ] **Step 1: Update Python regression script**

Modify `tests/regression_runner.py` to recursively iterate over all `.inp`/`.INP` files in `ISORROPIALite_Executable_Manual_Papers`, run both `isorropia_cli` and Fortran `isolite`, and assert relative accuracy.

- [ ] **Step 2: Verify entire integrated pipeline**

Run: `python3 tests/regression_runner.py`
Expected: 100% E2E validation PASS across all test configurations!

- [ ] **Step 3: Commit**

```bash
git add tests/regression_runner.py
git commit -m "test: scale E2E regression harness to check all forward and reverse input configurations"
```
