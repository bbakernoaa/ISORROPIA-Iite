include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(isorropia_tests test_main.cpp test_state.cpp test_thermo.cpp test_solver.cpp test_activities.cpp test_crustal.cpp test_reverse.cpp test_c_api.cpp test_properties.cpp test_smoothness.cpp test_batched.cpp)
target_link_libraries(isorropia_tests PRIVATE isorropia gtest_main)

include(GoogleTest)
gtest_discover_tests(isorropia_tests)

# Fortran Linkage Test
add_executable(isorropia_fortran_test test_linkage.f90)
target_link_libraries(isorropia_fortran_test PRIVATE isorropia)

add_test(NAME FortranLinkageTest COMMAND isorropia_fortran_test)
