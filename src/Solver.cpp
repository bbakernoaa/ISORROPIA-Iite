#include "Isorropia/Solver.hpp"

namespace Isorropia {

Solver::Solver() = default;

void Solver::solve(const Input& input, State& state) {
    // This is the main entry point (replaces 'SUBROUTINE ISOROPIA').
    // In Phase 1, we provide a placeholder that clears the error stack to verify the build.
    state.clear_errors();
}

} // namespace Isorropia
