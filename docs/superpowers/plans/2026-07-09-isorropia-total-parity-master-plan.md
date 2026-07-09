# ISORROPIA-Lite C++ Port (Absolute 1-to-1 Parity Master Plan)

This Master Plan outlines the systematic, phased roadmap to port **all remaining thermodynamic and numerical speciation subroutines** from the legacy Fortran codebase. This will eliminate all remaining skeletons and establish **absolute, 1-to-1 physical and chemical parity** across the entire ~20,000 lines of F77 code.

---

## 1. Context & Global Constraints

### The Porting Standard
- **No Skeletons**: Every single physical subroutine and iterative numerical step (Newton-Raphson solvers, dynamic root bisections, solid precipitate balances) must be completely translated to modernized, thread-safe C++17.
- **Double-Layered Verification**:
  1. **GoogleTest**: Unit-level assertions checking subroutines under isolated limits.
  2. **Python E2E Harness**: Side-by-side relative difference assertions comparing C++ outputs (`isorropia_cli`) directly to legacy Fortran outputs (`isolite`) on macOS:
     $$\max \left( \frac{|X_{C++} - X_{Fortran}|}{\max(|X_{Fortran}|, 1.0)} \right) \le 10^{-12}$$

---

## 2. Phased Roadmap to Total Parity

```
      ┌────────────────────────────────────────────────────────┐
      │  Phase 8a: Scaffolding, Constants, and Case 1 Solver   │  ◄── [COMPLETED]
      └───────────────────────────┬────────────────────────────┘
                                  │
                                  ▼
      ┌────────────────────────────────────────────────────────┐
      │  Phase 8b: Case 2 Speciation and Nitric Acid Sublim.   │
      └───────────────────────────┬────────────────────────────┘
                                  │
                                  ▼
      ┌────────────────────────────────────────────────────────┐
      │  Phase 8c: Crustal & Marine Speciation (Cases 3 & 4)   │
      └───────────────────────────┬────────────────────────────┘
                                  │
                                  ▼
      ┌────────────────────────────────────────────────────────┐
      │  Phase 8d: Standard Reverse Solvers (isrp1r - isrp4r)  │
      └───────────────────────────┬────────────────────────────┘
                                  │
                                  ▼
      ┌────────────────────────────────────────────────────────┐
      │  Phase 8e: Production Regression and Linking Approval  │
      └────────────────────────────────────────────────────────┘
```

---

### Phase 8b: Case 2 Speciation & Nitric Acid Sublimation

**Goal**: Port Speciation and sublimation equations for **Ammonium-Sulfate-Nitrate Metastable Aerosol systems (Case 2)** to achieve complete dynamic solvers coverage for all standard non-crustal configurations.

* **Task 1: Port Case D3 Speciation (`CALCD3` and `FUNCD3`)**
  * **Files**: `include/Isorropia/Solver.hpp`, `src/ForwardSolvers.cpp`
  * **Description**: Port the bisection search subroutine `CALCD3` and its root objective evaluator function `FUNCD3` (reproducing `isofwd.f` line 723). Iterates over hydrogen ion concentration ($H^+$) inside the liquid phase to solve coupled ammonium-nitrate equilibria.
* **Task 2: Port Gaseous Nitric Acid Dissolution Solver (`CALCNA`)**
  * **Files**: `include/Isorropia/Solver.hpp`, `src/ForwardSolvers.cpp`
  * **Description**: Port the dynamic sublimation equations for Nitric Acid (`CALCNA` in `isocom.f` line 3100), calculating gaseous partitions and back-adjusting H+ and nitrate ions in the liquid phase.
* **Task 3: Enable Strict E2E Regression Checks for `test1.inp`**
  * **Files**: `tests/regression_runner.py`
  * **Description**: Since all active chemistry for Cases 1 & 2 is now 100% active, promote `test1.inp` side-by-side relative difference assertions to strict double-precision bounds ($10^{-12}$).

---

### Phase 8c: Crustal & Marine Speciation (Cases 3 & 4)

**Goal**: Port the dynamic speciation and precipitation solvers for **Marine (NaCl-H2O)** and **Crustal (Na-Ca-K-Mg-Cl-H2O) systems** to achieve complete forward physical parity.

* **Task 1: Port Case 3 Speciation (`ISRP3F`)**
  * **Files**: `src/ForwardSolvers.cpp`, `include/Isorropia/Solver.hpp`
  * **Description**: Port the Newton-Raphson solvers and solid salt precipitation/crystallization algorithms for Marine systems (reproducing `isofwd.f` lines 800-1500).
* **Task 2: Port Case 4 Speciation (`ISRP4F`)**
  * **Files**: `src/ForwardSolvers.cpp`, `include/Isorropia/Solver.hpp`
  * **Description**: Port the massive, complex speciation, mineral-dissociation, and precipitation matrices for the full crustal configuration containing Calcium, Potassium, and Magnesium.
* **Task 3: Enable Strict E2E Checks for `Partitioning_with_organics.INP`**
  * **Files**: `tests/regression_runner.py`
  * **Description**: Promote the crustal E2E regression check to strict relative tolerance bounds ($10^{-12}$).

---

### Phase 8d: Standard Reverse Solvers (isrp1r - isrp4r)

**Goal**: Port the complete suite of **Reverse/Inverse Formulation Solvers** to achieve complete 1-to-1 reverse parity.

* **Task 1: Port Case 1 & 2 Reverse Solvers (`isrp1r`, `isrp2r`)**
  * **Files**: `src/ReverseSolvers.cpp`
  * **Description**: Port the dynamic speciation algorithms from `isorev.f` lines 1-1000 for standard systems in reverse mode.
* **Task 2: Port Case 3 & 4 Crustal/Marine Reverse Solvers (`isrp3r`, `isrp4r`)**
  * **Files**: `src/ReverseSolvers.cpp`
  * **Description**: Port the crustal, mineral, and marine reverse solvers from `isorev.f` lines 1000-1761.
* **Task 3: Enable Strict E2E Checks for `Reverse_with_organics.INP`**
  * **Files**: `tests/regression_runner.py`
  * **Description**: Promote reverse problem validations to strict relative tolerance bounds ($10^{-12}$).

---

### Phase 8e: Production Linking & CATChem Sign-off

**Goal**: Conduct exhaustive edge-case evaluations, optimize compilation flags, and export the library for production-grade CATChem integration.

* **Task 1: Standardize Linkage Exports**
  * **Files**: `CMakeLists.txt`
  * **Description**: Support standard CMake packaging (`IsorropiaConfig.cmake`), exporting library targets to simplify linking via `find_package(Isorropia REQUIRED)`.
* **Task 2: Compiler Profiling & Fast-Math Checks**
  * **Files**: `CMakeLists.txt`
  * **Description**: Profile compiler optimizations (`-O3`, `-march=native`) and check their influence on IEEE-754 numerical compliance to ensure absolute speed and numerical safety.
