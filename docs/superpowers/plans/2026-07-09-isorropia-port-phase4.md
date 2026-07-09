# ISORROPIA-Lite C++ Port (Phase 7: Modern Fortran 2003 Bindings Wrapper) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement a standard Fortran 2003 module using `ISO_C_BINDING` to map plain C structs and linkages to native Fortran type declarations, enabling seamless, high-performance integration with Fortran host models like CATChem.

**Architecture:** Create a Fortran 2003 module file `src/IsorropiaFortran.f90` that maps plain memory structures (`IsorropiaInput`, `IsorropiaState`) and link symbols via standard inter-language protocols, bypassing legacy compiler capitalizations.

**Tech Stack:** C++17, Fortran 2003 (gfortran), CMake

## Global Constraints

- Numerical Equivalence: Tolerance-based equality ($10^{-12}$ to $10^{-15}$).
- Thread Safety: No global variables or COMMON blocks.
- Robustness: Zero usage of C++ exceptions (`throw`) in solvers.
- Library Structure: Standalone C++ and Fortran compiled library using CMake.
- 0-Based Indexing for all C++ arrays.
- Scientist Readability: All variables, units, and chemical equations must be extensively documented in comments mapping directly back to original Fortran names and the scientific literature.

---

### Task 1: Fortran ISO_C_BINDING Wrapper Module

**Files:**
- Create: `src/IsorropiaFortran.f90`

**Interfaces:**
- Produces: Compiled Fortran 2003 module `isorropiafortran.mod` containing native handles and standard linkage interface bindings.

- [ ] **Step 1: Write src/IsorropiaFortran.f90**

```fortran
module IsorropiaFortran
  use, intrinsic :: iso_c_binding
  implicit none

  ! Modern Fortran types mapping exactly to standard Isorropia C layout structures
  type, bind(C) :: IsorropiaInput
    real(c_double) :: w(8)
    real(c_double) :: org(3)
    real(c_double) :: waer(8)
    real(c_double) :: temp
    real(c_double) :: rh
    integer(c_int) :: iprob
    integer(c_int) :: nadj
  end type IsorropiaInput

  type, bind(C) :: IsorropiaState
    real(c_double) :: temp
    real(c_double) :: rh
    real(c_double) :: w(8)
    real(c_double) :: waer(8)
    real(c_double) :: org(3)

    real(c_double) :: molal(10)
    real(c_double) :: molalr(23)
    real(c_double) :: gama(23)
    real(c_double) :: zz(23)
    real(c_double) :: z(10)
    real(c_double) :: gamou(23)
    real(c_double) :: gamin(23)
    real(c_double) :: m0(23)
    real(c_double) :: gasaq(3)
    integer(c_int) :: actmod
    real(c_double) :: epsact
    real(c_double) :: coh
    real(c_double) :: chno3
    real(c_double) :: chcl
    real(c_double) :: water
    real(c_double) :: ionic
    real(c_double) :: watcmp(24)
    integer(c_int) :: frst
    integer(c_int) :: calain
    integer(c_int) :: calaou
    integer(c_int) :: dryf

    real(c_double) :: ch2so4, cnh42s4, cnh4hs4, cnacl, cna2so4, cnano3, cnh4no3, cnh4cl, cnahso4, clc
    real(c_double) :: ccaso4, ccano32, ccacl2, ck2so4, ckhso4, ckno3, ckcl, cmgso4, cmgno32, cmgcl2
    real(c_double) :: gnh3, ghno3, ghcl

    integer(c_int) :: num_errors
  end type IsorropiaState

  ! Define standard interlanguage binding interfaces to link directly to flat C symbols
  interface
    subroutine isorropia_solve_c(input, state) bind(C, name="isorropia_solve_c")
      import :: IsorropiaInput, IsorropiaState
      type(IsorropiaInput), intent(in)  :: input
      type(IsorropiaState), intent(out) :: state
    end subroutine isorropia_solve_c
  end interface

contains

  ! Native high-level Fortran subroutine wrapping C solves for host integration (like CATChem)
  subroutine isorropia_solve_f(input, state)
    type(IsorropiaInput), intent(in)   :: input
    type(IsorropiaState), intent(out)  :: state
    call isorropia_solve_c(input, state)
  end subroutine isorropia_solve_f

end module IsorropiaFortran
```

- [ ] **Step 2: Commit**

```bash
git add src/IsorropiaFortran.f90
git commit -m "feat: implement modern Fortran 2003 wrapper module utilizing ISO_C_BINDING"
```

---

### Task 2: CMake Integration with Multi-language Compilation

**Files:**
- Modify: `CMakeLists.txt`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: C++ static targets and Fortran source.
- Produces: Standalone binary and link modules for tests.

- [ ] **Step 1: Enable Fortran compilation in root CMakeLists.txt**

Modify `project` declaration and targets:
```cmake
cmake_minimum_required(VERSION 3.14)
project(IsorropiaLite LANGUAGES CXX Fortran)

# ...
# Library target compiling both C++ and Fortran sources
add_library(isorropia STATIC 
  src/Solver.cpp 
  src/Thermodynamics.cpp 
  src/WaterActivities.cpp 
  src/ForwardSolvers.cpp 
  src/ActivityCoefficients_KM.cpp
  src/ReverseSolvers.cpp
  src/IsorropiaFortran.f90
)
```

- [ ] **Step 2: Commit**

```bash
git add CMakeLists.txt
git commit -m "chore: enable multi-language CMake compilation supporting CXX and Fortran"
```

---

### Task 3: Fortran Linkage Test Execution

**Files:**
- Create: `tests/test_linkage.f90`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `isorropia` module target.
- Produces: Compiled executable `isorropia_fortran_test` validating complete C++ to Fortran linking success.

- [ ] **Step 1: Write tests/test_linkage.f90**

```fortran
program test_linkage
  use IsorropiaFortran
  implicit none

  type(IsorropiaInput) :: input
  type(IsorropiaState) :: state

  print *, "=== Isorropia Fortran Linkage Test ==="

  ! Replicate Run 1 inputs
  input%w = 0.0d0
  input%w(2) = 1.0d0 ! H2SO4 component (1-indexed inside Fortran array)
  input%w(3) = 2.0d0 ! NH3 component
  input%w(4) = 1.0d0 ! HNO3 component

  input%org = 0.0d0
  input%org(1) = 10.0d0 ! org
  input%org(2) = 0.15d0 ! korg
  input%org(3) = 1000.0d0 ! density

  input%rh = 0.80d0
  input%temp = 298.15d0
  input%iprob = 0
  input%nadj = 1

  ! Call standard Fortran solver wrapper
  call isorropia_solve_f(input, state)

  print *, "Water Activity Content (kg/m3): ", state%water
  print *, "Diagnostic Errors:              ", state%num_errors

  if (state%water > 8.0d0 .and. state%num_errors == 0) then
    print *, "✅ FORTRAN-TO-C++ BINDING VERIFIED SUCCESSFULLY!"
  else
    print *, "❌ FORTRAN-TO-C++ BINDING LINKAGE TEST FAILED."
    call exit(1)
  end if

end program test_linkage
```

- [ ] **Step 2: Update tests/CMakeLists.txt to compile test binary**

Add to `tests/CMakeLists.txt`:
```cmake
# Fortran Linkage verification executable
add_executable(isorropia_fortran_test test_linkage.f90)
target_link_libraries(isorropia_fortran_test PRIVATE isorropia)

add_test(NAME FortranLinkageTest COMMAND isorropia_fortran_test)
```

- [ ] **Step 3: Build and run test suite**

Run: `cd build && cmake .. && make && ctest -V && python3 ../tests/regression_runner.py`
Expected: Passes `FortranLinkageTest` and all 19 C++ tests successfully!

- [ ] **Step 4: Commit**

```bash
git add tests/test_linkage.f90 tests/CMakeLists.txt
git commit -m "test: add Fortran linkage test checks compiled under CTest"
```
