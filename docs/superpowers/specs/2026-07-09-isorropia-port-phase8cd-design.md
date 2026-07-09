# ISORROPIA-Lite C++ Port: Phase 8c & 8d Speciation and Reverse Solvers Design

This document specifies the technical design, architectural routing, and structured loop conversions to port the remaining legacy Metastable (liquid-only) speciation solvers from Fortran 77 to C++17. This design achieves absolute 1-to-1 physical parity and removes all remaining skeletons from the C++ library.

---

## 1. Architectural Principles

### 1.1 Complete Absence of GOTO Statements
To comply with strict modern C++ guidelines, all legacy `GOTO` statements, multi-branching jumps, and interactive bisection branches are entirely banned. We will map them to structured flow-control mechanisms:
- **Root-Tracking Brackets:** Use structured scanning loops (`for` or `while`) over defined intervals to find sign-crossings.
- **Bisections:** Restructured as `while` loops with early `break` statements when convergence criteria (`std::abs(y2) <= eps`) are met, incrementing loop-bounds safely without back-jumps.
- **Exit Conditions:** Use Boolean flags (`bool success = false;`) and early `return` paths.

### 1.2 Thread Safety and Local State
To ensure absolute thread safety for multi-threaded or multi-cell integrations (like CATChem), all internal state is completely localized inside the `Isorropia::State` structure. No global parameters, static variables, or Fortran-like `COMMON` blocks are used.

### 1.3 0-Based Indexing and Variable Alignment
We strictly follow 0-based indexing for all arrays (e.g., mapping `MOLAL(1)` in Fortran to `molal[0]` in C++).

---

## 2. Phase 8c: Case 4 Crustal Forward Solvers (`isrp4f`)

The Crustal Forward solver manages Na-NH4-SO4-NO3-Cl-Ca-K-Mg-H2O systems in forward mode.

### 2.1 Driver Subroutine `Solver::isrp4f`
- **Replicates:** `SUBROUTINE ISRP4F` in `isofwd.f`.
- **Flow:**
  1. Initialize concentrations and compute dynamic deliquescence relative humidities.
  2. **Crustal Adjustments:** Check for excess cations. If `w[0]+w[6]+w[7]+w[8] > 2.0*w[1] + w[3] + w[4]`, sequentially calculate and subtract CaSO4, Ca(NO3)2, CaCl2, Na2SO4, NaCl, NaNO3, MgSO4, Mg(NO3)2, MgCl2, K2SO4, KHSO4, KNO3, KCl solid precipitations.
  3. **Ratio Evaluations:**
     - `so4rat = (w[0] + w[2] + w[5] + w[6] + w[7]) / w[1]`
     - `crnarat = (w[0] + w[5] + w[6] + w[7]) / w[1]`
     - `crrat = (w[5] + w[6] + w[7]) / w[1]`
  4. **Routing Dispatch:**
     - `so4rat >= 2.0 && crnarat < 2.0` $\rightarrow$ Case O7 (`cal_co7` and `funco7` evaluation)
     - `so4rat >= 2.0 && crnarat >= 2.0 && crrat <= 2.0` $\rightarrow$ Case M8 (`cal_cm8` and `funcm8` evaluation)
     - `so4rat >= 2.0 && crnarat >= 2.0 && crrat > 2.0` $\rightarrow$ Case P13 (`cal_cp13` and `funcp13` evaluation)
     - `1.0 <= so4rat < 2.0` $\rightarrow$ Case L9 (`cal_cl9` and dry `cal_cl1a`), followed by `cal_cnha` and `cal_cnh3`.
     - `so4rat < 1.0` $\rightarrow$ Case K4 (`cal_ck4` dry balance), followed by `cal_cnha` and `cal_cnh3`.

### 2.2 Case O7 (`cal_co7` and `funco7`)
- **Physics:** Sulfate poor, dust and sodium poor Metastable liquid.
- **Bisection:** Scan over the active bracket of $H^+$ ion concentrations. The objective function is electroneutrality.
- **Conversion:** Restructure the Fortran root-tracker (originally jumping back/forth) into:
  ```cpp
  double dx = (omehi - omelo) / static_cast<double>(ndiv);
  double x1 = omehi;
  double y1 = funco7(x1, input, state);
  double x2 = x1;
  double y2 = y1;
  // Step through intervals to find crossing
  for (int i = 1; i <= ndiv; ++i) {
      x2 = std::max(x1 - dx, omelo);
      y2 = funco7(x2, input, state);
      if (y1 * y2 < 0.0) {
          break; // sign change found
      }
      x1 = x2;
      y1 = y2;
  }
  // Bisection loop
  for (int it = 0; it < maxit; ++it) {
      double xmid = 0.5 * (x1 + x2);
      double ymid = funco7(xmid, input, state);
      if (std::abs(ymid) <= eps) {
          x2 = xmid;
          break;
      }
      if (y1 * ymid < 0.0) {
          x2 = xmid;
      } else {
          x1 = xmid;
          y1 = ymid;
      }
  }
  ```

### 2.3 Cases M8, P13, L9, K4
- **M8 & P13:** Run identical structured bisections over $H^+$ evaluating electroneutrality in dust-poor and dust-rich regimes respectively.
- **L9 & K4:** Completely analytical solvers that construct dry ionic matrices (e.g., CaSO4, MgSO4, Na2SO4, NH4HSO4, etc.) and calculate water uptake before calling gas-liquid partitioning.

---

## 3. Phase 8d: Standard Reverse Solvers (`isrp1r` - `isrp4r`)

Reverse solvers map ambient aerosol mass concentrations directly to solution properties.

### 3.1 Reverse Drivers
- `isrp1r` (NH4-SO4-H2O)
- `isrp2r` (NH4-SO4-NO3-H2O)
- `isrp3r` (Na-NH4-SO4-NO3-Cl-H2O)
- `isrp4r` (Na-NH4-SO4-NO3-Cl-Ca-K-Mg-H2O)

### 3.2 Calculation Paths
- **Wet Paths:** If relative humidity `rh` exceeds mutual deliquescence relative humidity, calculate limiting wet sulfate ratio (`SULRATW`). If the input ratio is poor (`SULRAT >= SULRATW`), dispatch to liquid-only bisection child solvers:
  - `isrp1r` $\rightarrow$ `cal_s2` (bisection on NH4+/H+ equilibrium)
  - `isrp2r` $\rightarrow$ `cal_n3` (bisection on HNO3/NH3/NH4+ equilibrium)
  - `isrp3r` $\rightarrow$ `cal_q5` / `cal_q1a` (marine bisections)
  - `isrp4r` $\rightarrow$ `cal_r6`, `cal_v7`, or `cal_u8` (crustal bisections depending on ratios)
- **Dry/Rich Paths:** If `SULRAT < SULRATW` (Sulfate Rich), run a corresponding Forward solver (e.g., `cal_cb4`, `cal_cc2`, `cal_ci6`, `cal_cj3`, `cal_cl9`, `cal_ck4`) on dry concentrations to determine physical aerosol state, and then back-calculate gas phase concentrations using `cal_cnhp` and `cal_cnh3p` partitioning.

### 3.3 Reverse Sub-solvers (`cal_s2`, `cal_n3`, `cal_q5`, `cal_r6`, `cal_v7`, `cal_u8`, `cal_w13`)
Each child reverse solver executes a bisection search inside liquid phase boundaries. Like Phase 8c, all back-jumping bisection logic will be restructured into clean, structured scan-and-bisection `while` loops.

---

## 4. Quality Control and Validation

### 4.1 Python End-to-End Regression Harness
All validations will run under `tests/regression_runner.py`. Once Phase 8c and 8d are implemented:
1. Promote `Partitioning_with_organics.INP` (Crustal E2E) from skeleton mode to strict comparison check.
2. Promote `Reverse_with_organics.INP` (Reverse E2E) from skeleton mode to strict comparison check.
3. Assert that relative differences are less than or equal to $10^{-12}$:
   $$\max \left( \frac{|X_{C++} - X_{Fortran}|}{\max(|X_{Fortran}|, 1.0)} \right) \le 10^{-12}$$

### 4.2 Unit Tests
Unit tests inside `tests/test_crustal.cpp` and `tests/test_reverse.cpp` will verify that appropriate solver routing occurs and that no internal diagnostics or convergence errors are flagged during runs.
