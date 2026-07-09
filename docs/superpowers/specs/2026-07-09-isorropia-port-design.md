# Design Spec: Porting Legacy ISORROPIA-Lite Fortran to Modernized C++

- **Date:** 2026-07-09
- **Status:** Proposed
- **Author:** Gemini CLI Agent

---

## 1. Context & Goals

**ISORROPIA-Lite** is a simplified and computationally optimized version of ISORROPIA-II, a widely-used atmospheric aerosol thermodynamics model that calculates the physical state (liquid/solid/gas) and chemical composition of aerosols.

The goal is to port the legacy Fortran codebase (~20,000 lines of code) to **modernized, thread-safe, and highly-performant C++** for integration into **CATChem**. 

### Critical Success Criteria
1. **Numerical Equivalence**: The C++ port must match the legacy Fortran results within a strict relative numerical tolerance of $10^{-12}$ to $10^{-15}$ (tolerance-based numerical equality). This accounts for minor differences in compiler math libraries, instruction reordering, and transcendental functions.
2. **Thread Safety**: Complete elimination of all global variables and `COMMON` blocks to allow thread-safe execution across CATChem's spatial grid cells.
3. **Robustness & Performance**: Zero usage of C++ exceptions (`throw`) in the simulation hot-path to avoid runtime overhead and prevent simulation crashes. A non-throwing, local error stack will be used instead.
4. **Verifiability**: Establish an automated, double-layered validation harness combining C++ Unit Tests (GoogleTest) and an end-to-end Python regression test suite comparing Fortran outputs to C++ outputs side-by-side.
5. **Scientist Readability**: High readability for atmospheric scientists. We must use clear, descriptive naming alongside comments mapping back to original Fortran variable names and peer-reviewed physical/chemical equations. All physical and chemical units must be explicitly documented in comments.

---

## 2. Architectural Decisions

### A. Library Structure
We will implement the port as a **Standalone C++ Compiled Library (C++17)**.
* **Why**: The codebase is large (~20,000 lines). A compiled library keeps compile times low for downstream consumers like CATChem, isolates implementation details, and simplifies unit test compilation.
* **Build System**: CMake.

### B. State Encapsulation
All legacy Fortran `COMMON` blocks will be fully encapsulated inside C++ data structures. 

```cpp
namespace Isorropia {

struct Input {
    std::array<double, 8> w = {0.0};       // W(NCOMP)
    std::array<double, 3> org = {0.0};     // ORG(NORG)
    std::array<double, 8> waer = {0.0};    // WAER(NCOMP)
    double temp = 298.15;                  // TEMP
    double rh = 0.0;                       // RH
    int iprob = 0;                         // IPROB
    int nadj = 0;                          // NADJ
};

struct ErrorEntry {
    int code = 0;
    std::string message;
};

struct State {
    // Replaces COMMON /IONS/
    std::array<double, 10> molal = {0.0};    // MOLAL(NIONS)
    std::array<double, 23> molalr = {0.0};   // MOLALR(NPAIR)
    std::array<double, 23> gama = {0.0};     // GAMA(NPAIR)
    std::array<double, 23> zz = {0.0};       // ZZ(NPAIR)
    std::array<double, 10> z = {0.0};        // Z(NIONS)
    std::array<double, 23> gamou = {0.0};    // GAMOU(NPAIR)
    std::array<double, 23> gamin = {0.0};    // GAMIN(NPAIR)
    std::array<double, 23> m0 = {0.0};       // M0(NPAIR)
    std::array<double, 3> gasaq = {0.0};     // GASAQ(NGASAQ)
    int actmod = 0;                          // ACTMOD
    double epsact = 0.0;
    double coh = 0.0;
    double chno3 = 0.0;
    double chcl = 0.0;
    double water = 0.0;
    float ionic = 0.0f;
    std::array<double, 24> watcmp = {0.0};   // WATCMP(NPAIR+1)
    bool frst = true;
    bool calain = false;
    bool calaou = false;
    bool dryf = false;

    // Replaces COMMON /SALT/
    double ch2so4 = 0.0, cnh42s4 = 0.0, cnh4hs4 = 0.0, cnacl = 0.0, cna2so4 = 0.0;
    double cnano3 = 0.0, cnh4no3 = 0.0, cnh4cl = 0.0, cnahso4 = 0.0, clc = 0.0;
    double ccaso4 = 0.0, ccano32 = 0.0, ccacl2 = 0.0, ck2so4 = 0.0, ckhso4 = 0.0;
    double ckno3 = 0.0, ckcl = 0.0, cmgso4 = 0.0, cmgno32 = 0.0, cmgcl2 = 0.0;

    // Replaces COMMON /GAS/
    double gnh3 = 0.0, ghno3 = 0.0, ghcl = 0.0;

    // Replaces COMMON /ZSR/ and /EQUK/
    std::array<double, 100> awas = {0.0}, awss = {0.0}, awac = {0.0}, awsc = {0.0};
    // ... all other ZSR activity arrays ...

    // Error Stack (Thread-Safe, Localized)
    bool stack_overflow = false;
    size_t num_errors = 0;
    std::array<ErrorEntry, 25> error_stack;

    void push_error(int code, std::string_view message) {
        if (num_errors >= error_stack.size()) {
            stack_overflow = true;
            return;
        }
        error_stack[num_errors] = {code, std::string(message)};
        num_errors++;
    }

    void clear_errors() {
        num_errors = 0;
        stack_overflow = false;
    }
};

class Solver {
public:
    Solver();
    void solve(const Input& input, State& state);
};

} // namespace Isorropia
```

---

## 3. Structural Organization & File Mapping

To make the ~20,000 lines of ported code maintainable, we map the original files to modern modular headers and source files:

| Source File | Content Category | Original Fortran Origin |
|---|---|---|
| `include/Isorropia/Solver.hpp` | Declarations of `Input`, `State`, `Solver` classes and enums. | `isrpia.inc` |
| `src/Solver.cpp` | Main entry point (`Solver::solve`), initialization logic, index mapping. | `main.f`, `isocom.f` |
| `src/Thermodynamics.cpp` | ZSR lookups, DRH calculations, chemical constants. | `isocom.f` |
| `src/ActivityCoefficients.cpp`| Pure salt activities, liquid activity coefficients solvers. | `isocom.f` |
| `src/ForwardSolvers.cpp` | Forward case sub-solvers (`ISRP1F`, `ISRP2F`, `ISRP3F`, `ISRP4F`). | `isofwd.f` |
| `src/ReverseSolvers.cpp` | Reverse case sub-solvers (`ISRP1R`, `ISRP2R`, `ISRP3R`, `ISRP4R`). | `isorev.f` |

---

## 4. Translation Specifics & Conventions

1. **0-Based Indexing Shifting**: 
   All arrays will be converted from Fortran 1-based indexing to C++ 0-based indexing.
2. **Descriptive Enums**:
   Enums will be defined to map array indices to chemical names:
   ```cpp
   enum class Component : size_t { Na = 0, H2SO4 = 1, NH3 = 2, HNO3 = 3, HCl = 4, Ca = 5, K = 6, Mg = 7 };
   enum class Ions : size_t { Na = 0, H = 1, NH4 = 2, NO3 = 3, Cl = 4, SO4 = 5, HSO4 = 6, Ca = 7, K = 8, Mg = 9 };
   ```
3. **Matrix Layouts**:
   Any 2D matrices will be converted from column-major to row-major format, with loop indices swapped accordingly to maintain CPU stride-1 cache performance.
4. **Error Handling**:
   No exceptions. Local error structures within the `State` class will log issues and allow the host application (CATChem) to handle non-converged grid cells gracefully.

---

## 5. Testing & Verification Plan

We will implement a dual-layer testing strategy:

```
                          ┌───────────────────────────┐
                          │     Porting Verification  │
                          └─────────────┬─────────────┘
                                        │
                ┌───────────────────────┴───────────────────────┐
                ▼                                               ▼
   ┌───────────────────────────┐                   ┌───────────────────────────┐
   │    C++ Unit Tests         │                   │   E2E Regression Tests    │
   │    (GoogleTest / GTest)   │                   │   (Python Script)         │
   ├───────────────────────────┤                   ├───────────────────────────┤
   │ * Pure salt activities    │                   │ * Runs Fortran & C++ bins │
   │ * DRH Calculations        │                   │ * Checks full input grid  │
   │ * State initialization    │                   │ * Asserts <= 1e-12 diff   │
   └───────────────────────────┘                   └───────────────────────────┘
```

1. **C++ Unit Tests (GTest)**:
   * Compiles under `tests/` directory.
   * Tests isolated components such as activity coefficients calculations, temperature corrections, and initial water activity lookups to guarantee correctness of internal equations.
2. **E2E Python Regression Suite**:
   * Uses input files like `test1.inp`.
   * Drives both the legacy Fortran binary and the newly compiled C++ executable, parses output files, and performs relative-difference verification.
   * Asserts:
     $$\max \left( \frac{|X_{C++} - X_{Fortran}|}{\max(|X_{Fortran}|, 1.0)} \right) \le 10^{-12}$$
