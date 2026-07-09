# ISORROPIA-Lite C++ Port (Phase 1: Infrastructure & Testing) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish the CMake build system, C++ core data structures (Input, State, Error stack), GoogleTest setup, and the E2E Python regression testing harness to prepare for the systematic porting of legacy Fortran solvers.

**Architecture:** A standalone C++17 library encapsulating legacy Fortran global state into a thread-safe `State` struct. Uses GTest for unit testing and a Python script for Fortran-vs-C++ output validation.

**Tech Stack:** C++17, CMake, GoogleTest, Python 3

## Global Constraints

- Numerical Equivalence: Tolerance-based equality ($10^{-12}$ to $10^{-15}$).
- Thread Safety: No global variables or COMMON blocks.
- Robustness: Zero usage of C++ exceptions (`throw`) in solvers.
- Library Structure: Standalone C++ Compiled Library (C++17) using CMake.
- 0-Based Indexing for all C++ arrays.
- Scientist Readability: All variables, units, and chemical equations must be extensively documented in comments mapping directly back to original Fortran names and the scientific literature.

---

### Task 1: Project Scaffolding and CMake Setup

**Files:**
- Create: `CMakeLists.txt`
- Create: `tests/CMakeLists.txt`
- Create: `tests/test_main.cpp`

**Interfaces:**
- Produces: CMake build system capable of building a static library `isorropia` and a test executable `isorropia_tests`.

- [ ] **Step 1: Write root CMakeLists.txt**

```cmake
cmake_minimum_required(VERSION 3.14)
project(IsorropiaLite CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library target
add_library(isorropia STATIC)
target_include_directories(isorropia PUBLIC include)

# Subdirectories
enable_testing()
add_subdirectory(tests)
```

- [ ] **Step 2: Write tests CMakeLists.txt with GTest**

```cmake
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(isorropia_tests test_main.cpp)
target_link_libraries(isorropia_tests PRIVATE isorropia gtest_main)

include(GoogleTest)
gtest_discover_tests(isorropia_tests)
```

- [ ] **Step 3: Write basic GTest main test**

```cpp
#include <gtest/gtest.h>

TEST(ScaffoldingTest, BasicAssertion) {
    EXPECT_TRUE(true);
}
```

- [ ] **Step 4: Verify build and test**

Run: `mkdir build && cd build && cmake .. && make && ctest --output-on-failure`
Expected: Passes `ScaffoldingTest.BasicAssertion`.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt tests/CMakeLists.txt tests/test_main.cpp
git commit -m "chore: scaffold CMake and GoogleTest"
```

---

### Task 2: Core Data Structures (Solver.hpp)

**Files:**
- Create: `include/Isorropia/Solver.hpp`
- Create: `src/Solver.cpp`

**Interfaces:**
- Consumes: Nothing
- Produces: `Isorropia::Input`, `Isorropia::State`, and `Isorropia::Solver` class stubs.

- [ ] **Step 1: Write Solver.hpp**

```cpp
#ifndef ISORROPIA_SOLVER_HPP
#define ISORROPIA_SOLVER_HPP

#include <array>
#include <string>
#include <string_view>

namespace Isorropia {

enum class Component : size_t { Na = 0, H2SO4 = 1, NH3 = 2, HNO3 = 3, HCl = 4, Ca = 5, K = 6, Mg = 7 };

struct Input {
    std::array<double, 8> w = {0.0};
    std::array<double, 3> org = {0.0};
    std::array<double, 8> waer = {0.0};
    double temp = 298.15;
    double rh = 0.0;
    int iprob = 0;
    int nadj = 0;
};

struct ErrorEntry {
    int code = 0;
    std::string message;
};

struct State {
    std::array<double, 10> molal = {0.0};
    std::array<double, 23> molalr = {0.0};
    std::array<double, 23> gama = {0.0};
    std::array<double, 23> zz = {0.0};
    std::array<double, 10> z = {0.0};
    std::array<double, 23> gamou = {0.0};
    std::array<double, 23> gamin = {0.0};
    std::array<double, 23> m0 = {0.0};
    std::array<double, 3> gasaq = {0.0};
    int actmod = 0;
    double epsact = 0.0;
    double coh = 0.0;
    double chno3 = 0.0;
    double chcl = 0.0;
    double water = 0.0;
    float ionic = 0.0f;
    std::array<double, 24> watcmp = {0.0};
    bool frst = true;
    bool calain = false;
    bool calaou = false;
    bool dryf = false;

    // Solid aerosol species
    double ch2so4 = 0.0, cnh42s4 = 0.0, cnh4hs4 = 0.0, cnacl = 0.0, cna2so4 = 0.0;
    double cnano3 = 0.0, cnh4no3 = 0.0, cnh4cl = 0.0, cnahso4 = 0.0, clc = 0.0;
    double ccaso4 = 0.0, ccano32 = 0.0, ccacl2 = 0.0, ck2so4 = 0.0, ckhso4 = 0.0;
    double ckno3 = 0.0, ckcl = 0.0, cmgso4 = 0.0, cmgno32 = 0.0, cmgcl2 = 0.0;

    // Gas species
    double gnh3 = 0.0, ghno3 = 0.0, ghcl = 0.0;

    // Error Stack
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

#endif
```

- [ ] **Step 2: Write Solver.cpp stub**

```cpp
#include "Isorropia/Solver.hpp"

namespace Isorropia {

Solver::Solver() = default;

void Solver::solve(const Input& input, State& state) {
    // Stub for future solver implementation
    state.clear_errors();
}

} // namespace Isorropia
```

- [ ] **Step 3: Update root CMakeLists.txt to include source**

```cmake
# Modify the existing target in CMakeLists.txt
target_sources(isorropia PRIVATE src/Solver.cpp)
```

- [ ] **Step 4: Verify build**

Run: `cd build && make`
Expected: Compiles successfully without errors.

- [ ] **Step 5: Commit**

```bash
git add include/Isorropia/Solver.hpp src/Solver.cpp CMakeLists.txt
git commit -m "feat: core state and input structures"
```

---

### Task 3: Unit Tests for State Error Stack

**Files:**
- Create: `tests/test_state.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Isorropia::State`

- [ ] **Step 1: Write test_state.cpp**

```cpp
#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"

TEST(StateTest, ErrorStackPushAndClear) {
    Isorropia::State state;
    
    state.push_error(1, "Test Error 1");
    state.push_error(2, "Test Error 2");
    
    EXPECT_EQ(state.num_errors, 2);
    EXPECT_FALSE(state.stack_overflow);
    EXPECT_EQ(state.error_stack[0].code, 1);
    EXPECT_EQ(state.error_stack[0].message, "Test Error 1");
    EXPECT_EQ(state.error_stack[1].code, 2);
    
    state.clear_errors();
    EXPECT_EQ(state.num_errors, 0);
    EXPECT_FALSE(state.stack_overflow);
}

TEST(StateTest, ErrorStackOverflow) {
    Isorropia::State state;
    
    for (int i = 0; i < 26; ++i) {
        state.push_error(i, "Overflow Test");
    }
    
    EXPECT_EQ(state.num_errors, 25);
    EXPECT_TRUE(state.stack_overflow);
}
```

- [ ] **Step 2: Update tests/CMakeLists.txt**

Modify `add_executable(isorropia_tests test_main.cpp)` to `add_executable(isorropia_tests test_main.cpp test_state.cpp)`

- [ ] **Step 3: Run tests**

Run: `cd build && make && ctest -V`
Expected: StateTest.ErrorStackPushAndClear and StateTest.ErrorStackOverflow PASS.

- [ ] **Step 4: Commit**

```bash
git add tests/test_state.cpp tests/CMakeLists.txt
git commit -m "test: add unit tests for state error stack"
```

---

### Task 4: Python E2E Regression Harness

**Files:**
- Create: `tests/regression_runner.py`

**Interfaces:**
- Consumes: Legacy Fortran outputs, C++ solver output (future task).

- [ ] **Step 1: Write Python script**

```python
import os
import sys

def parse_output_file(filepath):
    """Placeholder for parsing Fortran/C++ output files into a dict of values."""
    # In future, this will parse the specific values from ISORROPIA-Lite's text output
    return {"dummy_val": 1.0}

def compare_results(fortran_data, cpp_data, tolerance=1e-12):
    """Compares dictionaries of results."""
    for key in fortran_data:
        if key not in cpp_data:
            print(f"Key {key} missing in C++ output")
            return False
            
        f_val = fortran_data[key]
        c_val = cpp_data[key]
        
        diff = abs(c_val - f_val)
        rel_diff = diff / max(abs(f_val), 1.0)
        
        if rel_diff > tolerance:
            print(f"Mismatch for {key}: Fortran={f_val}, C++={c_val}, RelDiff={rel_diff}")
            return False
            
    return True

if __name__ == "__main__":
    print("Regression harness ready. Waiting for C++ binary implementation.")
    sys.exit(0)
```

- [ ] **Step 2: Test script execution**

Run: `python3 tests/regression_runner.py`
Expected: Output "Regression harness ready. Waiting for C++ binary implementation."

- [ ] **Step 3: Commit**

```bash
git add tests/regression_runner.py
git commit -m "test: setup python regression test harness"
```
