cmake_minimum_required(VERSION 3.14)
project(IsorropiaLite CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library target
add_library(isorropia STATIC src/Solver.cpp src/Thermodynamics.cpp src/WaterActivities.cpp src/ForwardSolvers.cpp)
target_include_directories(isorropia PUBLIC include)

# Executable target
add_executable(isorropia_cli src/main.cpp)
target_link_libraries(isorropia_cli PRIVATE isorropia)

# Subdirectories
enable_testing()
add_subdirectory(tests)
