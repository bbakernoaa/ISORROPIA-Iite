cmake_minimum_required(VERSION 3.14)
project(IsorropiaLite CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library target
add_library(isorropia STATIC src/Solver.cpp src/Thermodynamics.cpp)
target_include_directories(isorropia PUBLIC include)

# Subdirectories
enable_testing()
add_subdirectory(tests)
