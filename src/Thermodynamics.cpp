#include "Isorropia/Solver.hpp"
#include <cmath>
#include <algorithm>

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

void State::initialize_drh() {
    //-----------------------------------------------------------------------
    // Unicomponent DRH Defaults at 298.15K
    //-----------------------------------------------------------------------
    drh2so4  = 0.0;
    drnh42s4 = 0.7997;
    drnh4hs4 = 0.4000;
    drlc     = 0.6900;
    drnacl   = 0.7528;
    drnano3  = 0.7379;
    drnh4cl  = 0.7710;
    drnh4no3 = 0.6183;
    drna2so4 = 0.9300;
    drnahso4 = 0.5200;
    drcano32 = 0.4906;
    drcacl2  = 0.2830;
    drk2so4  = 0.9750;
    drkhso4  = 0.8600;
    drkno3   = 0.9248;
    drkcl    = 0.8426;
    drmgso4  = 0.8613;
    drmgno32 = 0.5400;
    drmgcl2  = 0.3284;

    //-----------------------------------------------------------------------
    // Temperature dependency corrections (if temperature is not 298.15K)
    // Replicates: IF (INT(TEMP) .NE. 298) in isocom.f
    //-----------------------------------------------------------------------
    if (static_cast<int>(temp) != 298) {
        double t0 = 298.15;
        double tcf = 1.0 / temp - 1.0 / t0;
        
        drnacl   *= std::exp(25.0 * tcf);
        drnano3  *= std::exp(304.0 * tcf);
        drna2so4 *= std::exp(80.0 * tcf);
        drnh4no3 *= std::exp(852.0 * tcf);
        drnh42s4 *= std::exp(80.0 * tcf);
        drnh4hs4 *= std::exp(384.0 * tcf);
        drlc     *= std::exp(186.0 * tcf);
        drnh4cl  *= std::exp(239.0 * tcf);
        drnahso4 *= std::exp(-45.0 * tcf);
        drcano32 *= std::exp(509.4 * tcf);
        drcacl2  *= std::exp(551.1 * tcf);
        drk2so4  *= std::exp(35.6 * tcf);
        drkcl    *= std::exp(159.0 * tcf);
        drmgso4  *= std::exp(-714.45 * tcf);
        drmgno32 *= std::exp(230.2 * tcf);
        drmgcl2  *= std::exp(42.23 * tcf);

        // Adjust for low-temperature DRH crossovers
        drnh4no3 = std::min({drnh4no3, drnh4cl, drnh42s4, drnano3, drnacl});
        drnano3  = std::min(drnano3, drnacl);
        drnh4cl  = std::min(drnh4cl, drnh42s4);
    }

    //-----------------------------------------------------------------------
    // Mutual Deliquescence Relative Humidities (MDRH) defaults
    //-----------------------------------------------------------------------
    drmlcab = 0.378;    // (NH4)3H(SO4)2 & NH4HSO4
    drmlcas = 0.690;    // (NH4)3H(SO4)2 & (NH4)2SO4
    drmasan = 0.600;    // (NH4)2SO4 & NH4NO3
    drmg1   = 0.460;    // (NH4)2SO4, NH4NO3, NA2SO4, NH4CL
    drmg2   = 0.691;    // (NH4)2SO4, NA2SO4, NH4CL
    drmg3   = 0.697;    // (NH4)2SO4, NA2SO4
    drmh1   = 0.240;    // NA2SO4, NANO3, NACL, NH4NO3, NH4CL
    drmh2   = 0.596;    // NA2SO4, NANO3, NACL, NH4CL
    drmi1   = 0.240;    // LC, NAHSO4, NH4HSO4, NA2SO4, (NH4)2SO4
    drmi2   = 0.363;    // LC, NAHSO4, NA2SO4, (NH4)2SO4
    drmi3   = 0.610;    // LC, NA2SO4, (NH4)2SO4
    drmq1   = 0.494;    // (NH4)2SO4, NH4NO3, NA2SO4
    drmr1   = 0.663;    // NA2SO4, NANO3, NACL
    drmr2   = 0.735;    // NA2SO4, NACL
    drmr3   = 0.673;    // NANO3, NACL
    drmr4   = 0.694;    // NA2SO4, NACL, NH4CL
    drmr5   = 0.731;    // NA2SO4, NH4CL
    drmr6   = 0.596;    // NA2SO4, NANO3, NH4CL
    drmr7   = 0.380;    // NA2SO4, NANO3, NACL, NH4NO3
    drmr8   = 0.380;    // NA2SO4, NACL, NH4NO3
    drmr9   = 0.494;    // NA2SO4, NH4NO3
    drmr10  = 0.476;    // NA2SO4, NANO3, NH4NO3
    drmr11  = 0.340;    // NA2SO4, NACL, NH4NO3, NH4CL
    drmr12  = 0.460;    // NA2SO4, NH4NO3, NH4CL
    drmr13  = 0.438;    // NA2SO4, NANO3, NH4NO3, NH4CL

    drmo1   = 0.460;    // (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4
    drmo2   = 0.691;    // (NH4)2SO4, NH4Cl, NA2SO4, K2SO4, MGSO4
    drmo3   = 0.697;    // (NH4)2SO4, NA2SO4, K2SO4, MGSO4
    drml1   = 0.240;    // K2SO4, MGSO4, KHSO4, NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
    drml2   = 0.363;    // K2SO4, MGSO4, KHSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
    drml3   = 0.610;    // K2SO4, MGSO4, KHSO4, (NH4)2SO4, NA2SO4, LC
    drmmm1  = 0.240;    // K2SO4, NA2SO4, MGSO4, NH4CL, NACL, NANO3, NH4NO3
    drmmm2  = 0.596;    // K2SO4, NA2SO4, MGSO4, NH4CL, NACL, NANO3
    drmp1   = 0.200;    // CA(NO3)2, CACL2, K2SO4, KNO3, KCL, ...
    drmp2   = 0.240;    // CA(NO3)2, K2SO4, KNO3, KCL, ...
    drmp3   = 0.240;    // CA(NO3)2, K2SO4, KNO3, KCL, ...
    drmp4   = 0.240;    // K2SO4, KNO3, KCL, ...
    drmp5   = 0.240;    // K2SO4, KNO3, KCL, ...
    drmv1   = 0.494;    // (NH4)2SO4, NH4NO3, NA2SO4, K2SO4, MGSO4
}

} // namespace Isorropia
