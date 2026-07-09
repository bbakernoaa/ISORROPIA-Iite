# ISORROPIA-Lite C++ Port: Phase 8c & 8d Speciation and Reverse Solvers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port Case 4 Forward speciation solvers and all Reverse solvers (isrp1r-isrp4r) from F77 to C++17, ensuring absolute 1-to-1 physical parity with zero `goto` statements.

**Architecture:** Extend `src/ForwardSolvers.cpp` and `src/ReverseSolvers.cpp` with the respective child sub-solvers. All state properties are completely localized inside the thread-safe `State` object. Restructure bisection searches using structured flow-control (`while` and `for` loops, early `break` statements).

**Tech Stack:** C++17, GoogleTest, CMake, Python 3

## Global Constraints
- Absolute numerical parity ($10^{-12}$ relative difference tolerance) compared to F77 on all outputs.
- No `goto` statements are allowed. All logic must use structured loop and branching flow-control.
- Strict 0-based array indexing for all mapped variables.
- All internal simulation state must be thread-safe (stored strictly inside `Isorropia::State`).
- Extensive inline commenting mapping variable names directly back to F77 and scientific documentation.

---

### Task 1: Declare all Case 4 Forward and Reverse Solvers in headers

**Files:**
- Modify: `include/Isorropia/Solver.hpp`

**Interfaces:**
- Declares private helper methods for Case 4 Forward speciation and all Reverse child solvers.

- [ ] **Step 1: Declare helpers in include/Isorropia/Solver.hpp**

Modify `include/Isorropia/Solver.hpp` to declare the following private member functions of `Solver` class:
```cpp
    // Case 4 Forward Solvers
    void cal_co7(const Input& input, State& state);
    double funco7(double x, const Input& input, State& state);
    void cal_cm8(const Input& input, State& state);
    double funcm8(double x, const Input& input, State& state);
    void cal_cp13(const Input& input, State& state);
    double funcp13(double x, const Input& input, State& state);
    void cal_cl9(const Input& input, State& state);
    void cal_cl1a(const Input& input, State& state);
    void cal_ck4(const Input& input, State& state);

    // Reverse Child Solvers
    void cal_s2(const Input& input, State& state);
    double funcs2(double x, const Input& input, State& state);
    void cal_n3(const Input& input, State& state);
    double funcn3(double x, const Input& input, State& state);
    void cal_q5(const Input& input, State& state);
    double funcq5(double x, const Input& input, State& state);
    void cal_q1a(const Input& input, State& state);
    void cal_r6(const Input& input, State& state);
    double funcr6(double x, const Input& input, State& state);
    void cal_r1a(const Input& input, State& state);
    void cal_v7(const Input& input, State& state);
    double funcv7(double x, const Input& input, State& state);
    void cal_v1a(const Input& input, State& state);
    void cal_u8(const Input& input, State& state);
    double funcu8(double x, const Input& input, State& state);
    void cal_u1a(const Input& input, State& state);
    void cal_w13(const Input& input, State& state);
    double funcw13(double x, const Input& input, State& state);
    void cal_w1a(const Input& input, State& state);
    void cal_nap(State& state);
```

- [ ] **Step 2: Re-compile project to ensure header updates compile successfully**

Run: `cmake -B build -S . && cmake --build build`
Expected: Success

- [ ] **Step 3: Commit declarations**

```bash
git add include/Isorropia/Solver.hpp
git commit -m "feat: declare Case 4 Forward and Reverse solver subroutines in Solver.hpp"
```

---

### Task 2: Implement Case 4 Forward Solvers (`isrp4f`, Case O7)

**Files:**
- Modify: `src/ForwardSolvers.cpp`

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: `void Solver::isrp4f(const Input&, State&)` and `void Solver::cal_co7(const Input&, State&)`

- [ ] **Step 1: Implement Excess Cations adjustment in Solver::isrp4f**

Port the sequential precipitation checks from `isofwd.f` line 265 into `isrp4f` to adjust cations (Ca, Na, Mg, K) when in excess. Then calculate the ratios `so4rat`, `crnarat`, `crrat` and dispatch to O7, M8, P13, L9, K4 solvers.
```cpp
void Solver::isrp4f(const Input& input, State& state) {
    state.clear_errors();
    state.scase = "4F";
    state.actmod = 4;

    // Port sequential precipitation adjustments here ...
    // ...
    double rest = 2.0 * state.w[1] + state.w[3] + state.w[4];
    if (state.w[0] + state.w[5] + state.w[6] + state.w[7] > rest) {
        // ... adjust as per isofwd.f ...
    }

    double so4rat  = (state.w[0] + state.w[2] + state.w[5] + state.w[6] + state.w[7]) / state.w[1];
    double crnarat = (state.w[0] + state.w[5] + state.w[6] + state.w[7]) / state.w[1];
    double crrat   = (state.w[5] + state.w[6] + state.w[7]) / state.w[1];

    if (so4rat >= 2.0 && crnarat < 2.0) {
        state.scase = "O7";
        cal_co7(input, state);
    } else if (so4rat >= 2.0 && crnarat >= 2.0 && crrat <= 2.0) {
        state.scase = "M8";
        cal_cm8(input, state);
    } else if (so4rat >= 2.0 && crnarat >= 2.0 && crrat > 2.0) {
        state.scase = "P13";
        cal_cp13(input, state);
    } else if (so4rat >= 1.0 && so4rat < 2.0) {
        state.scase = "L9";
        cal_cl9(input, state);
        cal_cnha(input, state);
        cal_cnh3(input, state);
    } else {
        state.scase = "K4";
        cal_ck4(input, state);
        cal_cnha(input, state);
        cal_cnh3(input, state);
    }
}
```

- [ ] **Step 2: Implement Case O7 Speciation (`cal_co7` and `funco7`)**

Replicate bisection search over Chloride evaporation `psi6` without `goto` statements.
```cpp
void Solver::cal_co7(const Input& input, State& state) {
    state.calaou = true;
    // Port CALCO7 ...
    // Restructure interval scan and bisection as detailed in Design spec Section 2.2
}
double Solver::funco7(double x, const Input& input, State& state) {
    // Port FUNCO7 replicating isofwd.f line 1790 ...
    // Uses state.cal_cmr() and state.cal_act2()
}
```

- [ ] **Step 3: Compile and run test_crustal to verify O7 routing**

Run: `cmake --build build && ./build/tests/test_crustal`
Expected: Runs successfully and routes appropriately.

- [ ] **Step 4: Commit Case O7 implementation**

```bash
git add src/ForwardSolvers.cpp
git commit -m "feat: implement isrp4f driver and Case O7 Forward speciation"
```

---

### Task 3: Implement Cases M8 and P13 Crustal Forward Solvers

**Files:**
- Modify: `src/ForwardSolvers.cpp`

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: `cal_cm8`, `funcm8`, `cal_cp13`, `funcp13`

- [ ] **Step 1: Port Case M8 speciation and objective function**

Translate `CALCM8` and `FUNCM8` from `isofwd.f` line 1922 and line 2024 respectively. Restructure bisection searches over $H^+$ ion concentration `psi6` as structured loops.

- [ ] **Step 2: Port Case P13 speciation and objective function**

Translate `CALCP13` and `FUNCP13` from `isofwd.f` line 2161 and line 2286 respectively. Restructure bisection searches as structured loops.

- [ ] **Step 3: Compile and run unit tests**

Run: `cmake --build build && ./build/tests/test_crustal`
Expected: Success

- [ ] **Step 4: Commit M8 and P13 solvers**

```bash
git add src/ForwardSolvers.cpp
git commit -m "feat: implement Case M8 and P13 Forward bisection solvers"
```

---

### Task 4: Implement Case L9 and K4 Crustal Forward Solvers

**Files:**
- Modify: `src/ForwardSolvers.cpp`

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: `cal_cl9`, `cal_cl1a`, `cal_ck4`

- [ ] **Step 1: Port Case L9 Speciation and dry material balances**

Translate `CALCL9` and `CALCL1A` from `isofwd.f` line 2459 and line 2564 respectively. This is an analytical solver setting solid matrices (CaSO4, CK2SO4, etc.) and liquid concentrations before doing ZSR water content updates.

- [ ] **Step 2: Port Case K4 Speciation**

Translate `CALCK4` from `isofwd.f` line 2651.

- [ ] **Step 3: Compile and verify execution**

Run: `cmake --build build && ./build/tests/test_crustal`
Expected: Success

- [ ] **Step 4: Commit L9 and K4 solvers**

```bash
git add src/ForwardSolvers.cpp
git commit -m "feat: implement Case L9 and K4 analytical Forward solvers"
```

---

### Task 5: Implement Reverse Drivers and Cases S2 & N3

**Files:**
- Modify: `src/ReverseSolvers.cpp`

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: Reverse drivers (`isrp1r`–`isrp4r`), `cal_s2`, `funcs2`, `cal_n3`, `funcn3`, `cal_nap`

- [ ] **Step 1: Port Case S2 reverse speciation and bisection**

Translate `CALCS2` from `isorev.f` line 441. Implement `funcs2` electroneutrality objective function. Restructure bisection search into a structured loop.

- [ ] **Step 2: Port Case N3 reverse speciation and bisection**

Translate `CALCN3` from `isorev.f` line 525. Implement `funcn3` objective function.

- [ ] **Step 3: Port Inverse Gas Nitric dissolution partitioner `CALCNAP`**

Translate `CALCNAP` from `isocom.f` line 2974 into `Solver::cal_nap(State& state)`. Sets up `ghno3` and active aquated species `gasaq[2]`.

- [ ] **Step 4: Update isrp1r and isrp2r reverse drivers to execute true physical solvers**

Modify `isrp1r` and `isrp2r` in `src/ReverseSolvers.cpp` to calculate `SULRATW` and conditionally route to `cal_s2` / `cal_n3` (wet paths) or forward solvers + reverse gas partitioning (dry/rich paths) as described in Design spec Section 3.2.

- [ ] **Step 5: Compile and run test_reverse to verify routing**

Run: `cmake --build build && ./build/tests/test_reverse`
Expected: Success

- [ ] **Step 6: Commit S2, N3, and drivers**

```bash
git add src/ReverseSolvers.cpp
git commit -m "feat: implement reverse drivers isrp1r-isrp2r and Case S2/N3 reverse solvers"
```

---

### Task 6: Implement Case Q5/Q1A Marine Reverse Solvers

**Files:**
- Modify: `src/ReverseSolvers.cpp`

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: `cal_q5`, `funcq5`, `cal_q1a`

- [ ] **Step 1: Port Case Q5 and Q1A speciation**

Translate `CALCQ5` and `CALCQ1A` from `isorev.f` line 639 and line 790 respectively. Restructure the bisection search over $H^+$ concentration as a structured loop.

- [ ] **Step 2: Update isrp3r reverse driver**

Update `isrp3r` in `src/ReverseSolvers.cpp` to calculate `SULRATW` and route to `cal_q5` or forward `isrp3f` + partitioners.

- [ ] **Step 3: Compile and run tests**

Run: `cmake --build build && ./build/tests/test_reverse`
Expected: Success

- [ ] **Step 4: Commit Case Q5/Q1A solvers**

```bash
git add src/ReverseSolvers.cpp
git commit -m "feat: implement Case Q5 and Q1A Marine Reverse solvers and update isrp3r driver"
```

---

### Task 7: Implement Crustal Reverse Solvers (Cases R6, V7, U8, W13)

**Files:**
- Modify: `src/ReverseSolvers.cpp`

**Interfaces:**
- Consumes: `Isorropia::State`
- Produces: `cal_r6`, `funcr6`, `cal_r1a`, `cal_v7`, `funcv7`, `cal_v1a`, `cal_u8`, `funcu8`, `cal_u1a`, `cal_w13`, `funcw13`, `cal_w1a`

- [ ] **Step 1: Port Case R6 & R1A speciation**

Translate `CALCR6` and `CALCR1A` from `isorev.f` line 838 and line 989 respectively.

- [ ] **Step 2: Port Case V7 & V1A speciation**

Translate `CALCV7` and `CALCV1A` from `isorev.f` line 1043 and line 1218 respectively.

- [ ] **Step 3: Port Case U8 & U1A speciation**

Translate `CALCU8` and `CALCU1A` from `isorev.f` line 1275 and line 1455 respectively.

- [ ] **Step 4: Port Case W13 & W1A speciation**

Translate `CALCW13` and `CALCW1A` from `isorev.f` line 1523 and line 1711 respectively.

- [ ] **Step 5: Update isrp4r reverse driver**

Update `isrp4r` in `src/ReverseSolvers.cpp` to calculate ratios and dispatch to `cal_r6`, `cal_v7`, `cal_u8`, `cal_w13` or dry forward solvers.

- [ ] **Step 6: Compile and run unit tests**

Run: `cmake --build build && ./build/tests/test_reverse`
Expected: Success

- [ ] **Step 7: Commit Crustal Reverse solvers**

```bash
git add src/ReverseSolvers.cpp
git commit -m "feat: implement Crustal Reverse solvers R6, V7, U8, W13 and update isrp4r driver"
```

---

### Task 8: Enable strict E2E comparisons and verify parity

**Files:**
- Modify: `tests/regression_runner.py`

**Interfaces:**
- Enables relative difference comparisons for Crustal and Reverse E2E configurations.

- [ ] **Step 1: Update tests/regression_runner.py to execute strict comparisons**

Modify `tests/regression_runner.py` line 202 to remove skeleton skips and enforce `compare_results` with relative tolerance of `1e-11` (or `1e-12`) for all input files: `test1.inp`, `Partitioning_with_organics.INP`, and `Reverse_with_organics.INP`.

- [ ] **Step 2: Execute build and run the regression runner**

Run: `cmake --build build && python tests/regression_runner.py`
Expected: ALL SELECTION INPUT SIMULATIONS REGRESSION VALIDATED SUCCESSFULLY! with PASS on all three files.

- [ ] **Step 3: Run full GoogleTest suite to ensure no regressions**

Run: `ctest --test-dir build`
Expected: 100% tests passed.

- [ ] **Step 4: Commit E2E test promotions**

```bash
git add tests/regression_runner.py
git commit -m "test: promote Partitioning_with_organics and Reverse_with_organics to strict E2E validation checks"
```
