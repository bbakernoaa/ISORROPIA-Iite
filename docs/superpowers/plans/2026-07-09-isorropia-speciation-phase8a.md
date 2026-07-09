# ISORROPIA-Lite C++ Port (Phase 8: True Case 1 Metastable Speciation Solver) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the physical speciation equations and iterative numerical solvers for pure NH4-SO4-H2O Metastable systems (originally subroutines `CALCA2`, `CALCB4`, `CALCC2`, and `CALCNH3` in `isofwd.f` and `isocom.f`), replacing the skeleton with live thermodynamics.

**Architecture:** Integrate bisection root finding and analytical speciation steps inside `src/ForwardSolvers.cpp`. These use the thread-local state and local error logging rather than exceptions.

**Tech Stack:** C++17, GoogleTest, CMake

## Global Constraints

- Numerical Equivalence: Tolerance-based equality ($10^{-12}$ to $10^{-15}$).
- Thread Safety: No global variables or COMMON blocks.
- Robustness: Zero usage of C++ exceptions (`throw`) in solvers.
- Library Structure: Standalone C++ Compiled Library (C++17) using CMake.
- 0-Based Indexing for all C++ arrays.
- Scientist Readability: All variables, units, and chemical equations must be extensively documented in comments mapping directly back to original Fortran names and the scientific literature.

---

### Task 1: Port Case A2 Speciation (Sulfate Poor Metastable)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare helper methods)
- Modify: `src/ForwardSolvers.cpp` (implement cal_ca2 and funca2 bisection equations)

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: Speciation calculations for Case A2.

- [ ] **Step 1: Declare helper functions inside Solver.hpp**

```cpp
private:
    void cal_ca2(const Input& input, State& state);
    double funca2(double x, const Input& input, State& state);
```

- [ ] **Step 2: Port CALCA2 and FUNCA2 from isofwd.f**

Replicate the bisection solver loop from `isofwd.f` line 400. Perform up to 100 bisection loops to find the roots of `funca2` within accuracy limits `eps` ($10^{-6}$).

---

### Task 2: Port Case B4 and C2 Speciation (Sulfate Rich)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare cal_cb4, cal_cc2)
- Modify: `src/ForwardSolvers.cpp` (implement speciation equations)

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: `void Solver::cal_cb4(const Input&, State&)` and `void Solver::cal_cc2(const Input&, State&)`

- [ ] **Step 1: Declare subroutines inside Solver.hpp**

```cpp
private:
    void cal_cb4(const Input& input, State& state);
    void cal_cc2(const Input& input, State& state);
```

- [ ] **Step 2: Implement cal_cb4 and cal_cc2 in src/ForwardSolvers.cpp**

Port the analytical speciation math for sulfate rich regimes from `isofwd.f` line 548 and line 668.

---

### Task 3: Port Gas-Liquid Ammonia Solver (CALCNH3)

**Files:**
- Modify: `include/Isorropia/Solver.hpp` (declare cal_cnh3)
- Modify: `src/ForwardSolvers.cpp` (implement cal_cnh3)

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: `void Solver::cal_cnh3(const Input&, State&)` setting `state.gnh3` and `state.gasaq`.

- [ ] **Step 1: Declare and implement cal_cnh3**

Port calculations from `isocom.f` line 3001, evaluating gaseous-liquid partition ratios.

- [ ] **Step 2: Replace isrp1f skeleton with live solver calls**

In `src/ForwardSolvers.cpp`, update `isrp1f` to execute the true physical tree:
```cpp
void Solver::isrp1f(const Input& input, State& state) {
    state.clear_errors();
    double sulrat = input.w[2] / input.w[1];
    
    if (sulrat >= 2.0) {
        state.scase = "A2";
        cal_ca2(input, state);
    } else if (sulrat >= 1.0) {
        state.scase = "B4";
        cal_cb4(input, state);
        cal_cnh3(input, state);
    } else {
        state.scase = "C2";
        cal_cc2(input, state);
        cal_cnh3(input, state);
    }
}
```

---

### Task 4: Complete Unit Testing and Verification

**Files:**
- Modify: `tests/test_solver.cpp`

- [ ] **Step 1: Run complete build and verify Case 1 Metastables pass**

Run: `cd build && make && ctest -V && python3 ../tests/regression_runner.py`
Expected: 100% Success!
