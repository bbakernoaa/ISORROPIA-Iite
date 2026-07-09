cmake_minimum_required(VERSION 3.14)
project(IsorropiaLite CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Generate a dummy source file in the build directory to satisfy CMake's requirement
# for at least one source file in a static library before target_sources is used in Task 2.
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp "")

# Library target
add_library(isorropia STATIC ${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp)
target_include_directories(isorropia PUBLIC include)

# Subdirectories
enable_testing()
add_subdirectory(tests)
