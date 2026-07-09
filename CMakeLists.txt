cmake_minimum_required(VERSION 3.14)
project(IsorropiaLite LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library target (C++ and Fortran sources compiled together)
add_library(isorropia STATIC 
  src/Solver.cpp 
  src/Thermodynamics.cpp 
  src/WaterActivities.cpp 
  src/ForwardSolvers.cpp 
  src/ActivityCoefficients_KM.cpp
  src/ReverseSolvers.cpp
  src/IsorropiaFortran.f90
)
target_include_directories(isorropia PUBLIC include)

# Executable target
add_executable(isorropia_cli src/main.cpp)
target_link_libraries(isorropia_cli PRIVATE isorropia)

# Subdirectories
enable_testing()
add_subdirectory(tests)
