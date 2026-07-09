#include "Isorropia/Solver.hpp"

namespace Isorropia {

void State::initialize_constants() {
    r = 82.0567e-6;
    tiny = 1.0e-20;
    tiny2 = 1.0e-11;
    great = 1.0e10;
    zero = 0.0;
    one = 1.0;

    // IMW initialization (NIONS=10)
    // Na+, H+, NH4+, NO3-, Cl-, SO4--, HSO4-, Ca++, K+, Mg++
    imw = {23.0, 1.0, 18.0, 62.0, 35.5, 96.0, 97.0, 40.1, 39.1, 24.3};

    // WMW initialization (NCOMP=8)
    // Na, H2SO4, NH3, HNO3, HCl, Ca, K, Mg
    wmw = {23.0, 98.0, 17.0, 63.0, 36.5, 40.1, 39.1, 24.3};

    // SMW initialization (NPAIR=23)
    smw = {
        85.0,  80.0,  58.5,  53.5,  142.0, 132.0, 120.0, 115.0, 120.0, 136.0,
        164.0, 111.0, 174.0, 136.0, 101.0, 74.5,  120.0, 148.0, 95.0,  18.0,
        0.0,   0.0,   0.0
    };
}

} // namespace Isorropia
