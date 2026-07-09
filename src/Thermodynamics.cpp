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

void State::calculate_equilibrium_constants() {
    //-----------------------------------------------------------------------
    // Base equilibrium constants at T0 = 298.15 K
    // Replicates INIT4 assignments in isocom.f
    //-----------------------------------------------------------------------
    xk1  = 1.015e-2;  // HSO4(aq)         <==> H(aq)     + SO4(aq)
    xk21 = 57.639;    // NH3(g)           <==> NH3(aq)
    xk22 = 1.805e-5;  // NH3(aq)          <==> NH4(aq)   + OH(aq)
    xk3  = 1.971e6;   // HCL(g)           <==> H(aq)     + CL(aq)
    xk31 = 2.500e3;   // HCL(g)           <==> HCL(aq)
    xk4  = 2.511e6;   // HNO3(g)          <==> H(aq)     + NO3(aq)
    xk41 = 2.100e5;   // HNO3(g)          <==> HNO3(aq)
    xk5  = 0.4799;    // NA2SO4(s)        <==> 2*NA(aq)  + SO4(aq)
    xk6  = 1.086e-16; // NH4CL(s)         <==> NH3(g)    + HCL(g)
    xk7  = 1.817;     // (NH4)2SO4(s)     <==> 2*NH4(aq) + SO4(aq)
    xk8  = 37.661;    // NACL(s)          <==> NA(aq)    + CL(aq)
    xk10 = 4.199e-17; // NH4NO3(s)        <==> NH3(g)    + HNO3(g) (Mozurkewich, 1993)
    xk11 = 2.413e4;   // NAHSO4(s)        <==> NA(aq)    + HSO4(aq)
    xk12 = 1.382e2;   // NH4HSO4(s)       <==> NH4(aq)   + HSO4(aq)
    xk13 = 29.268;    // (NH4)3H(SO4)2(s) <==> 3*NH4(aq) + HSO4(aq) + SO4(aq)
    xk14 = 22.05;     // NH4CL(s)         <==> NH4(aq)   + CL(aq)
    xkw  = 1.010e-14; // H2O              <==> H(aq)     + OH(aq)
    xk9  = 11.977;    // NANO3(s)         <==> NA(aq)    + NO3(aq)
    
    xk15 = 6.067e5;   // CA(NO3)2(s)      <==> CA(aq)    + 2NO3(aq)
    xk16 = 7.974e11;  // CACL2(s)         <==> CA(aq)    + 2CL(aq)
    xk17 = 1.569e-2;  // K2SO4(s)         <==> 2K(aq)    + SO4(aq)
    xk18 = 24.016;    // KHSO4(s)         <==> K(aq)     + HSO4(aq)
    xk19 = 0.872;     // KNO3(s)          <==> K(aq)     + NO3(aq)
    xk20 = 8.680;     // KCL(s)           <==> K(aq)     + CL(aq)
    xk23 = 1.079e5;   // MGS04(s)         <==> MG(aq)    + SO4(aq)
    xk24 = 2.507e15;  // MG(NO3)2(s)      <==> MG(aq)    + 2NO3(aq)
    xk25 = 9.557e21;  // MGCL2(s)         <==> MG(aq)    + 2CL(aq)

    //-----------------------------------------------------------------------
    // Temperature corrections (if temperature is not 298.15K)
    // Replicates van 't Hoff corrections in isocom.f
    //-----------------------------------------------------------------------
    if (static_cast<int>(temp) != 298) {
        double t0  = 298.15;
        double t0t = t0 / temp;
        double coef = 1.0 + std::log(t0t) - t0t;

        xk1  *= std::exp(8.85 * (t0t - 1.0) + 25.140 * coef);
        xk21 *= std::exp(13.79 * (t0t - 1.0) - 5.393 * coef);
        xk22 *= std::exp(-1.50 * (t0t - 1.0) + 26.920 * coef);
        xk3  *= std::exp(30.20 * (t0t - 1.0) + 19.910 * coef);
        xk31 *= std::exp(30.20 * (t0t - 1.0) + 19.910 * coef);
        xk4  *= std::exp(29.17 * (t0t - 1.0) + 16.830 * coef);
        xk41 *= std::exp(29.17 * (t0t - 1.0) + 16.830 * coef);
        xk5  *= std::exp(0.98 * (t0t - 1.0) + 39.500 * coef);
        xk6  *= std::exp(-71.00 * (t0t - 1.0) + 2.400 * coef);
        xk7  *= std::exp(-2.65 * (t0t - 1.0) + 38.570 * coef);
        xk8  *= std::exp(-1.56 * (t0t - 1.0) + 16.900 * coef);
        xk9  *= std::exp(-8.22 * (t0t - 1.0) + 16.010 * coef);
        xk10 *= std::exp(-74.7351 * (t0t - 1.0) + 6.025 * coef);
        xk11 *= std::exp(0.79 * (t0t - 1.0) + 14.746 * coef);
        xk12 *= std::exp(-2.87 * (t0t - 1.0) + 15.830 * coef);
        xk13 *= std::exp(-5.19 * (t0t - 1.0) + 54.400 * coef);
        xk14 *= std::exp(24.55 * (t0t - 1.0) + 16.900 * coef);
        xkw  *= std::exp(-22.52 * (t0t - 1.0) + 26.920 * coef);
        
        xk17 *= std::exp(-9.585 * (t0t - 1.0) + 45.81 * coef);
        xk18 *= std::exp(-8.423 * (t0t - 1.0) + 17.96 * coef);
        xk19 *= std::exp(-14.08 * (t0t - 1.0) + 19.39 * coef);
        xk20 *= std::exp(-6.902 * (t0t - 1.0) + 19.95 * coef);
        
        // Note: xk15, xk16, xk23, xk24, xk25 temperature corrections are 0.0 in Fortran
    }

    // Derived equilibrium constants
    xk2  = xk21 * xk22;
    xk42 = xk4 / xk41;
    xk32 = xk3 / xk31;
}

} // namespace Isorropia
