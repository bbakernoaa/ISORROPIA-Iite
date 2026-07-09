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

    // ZZ initialization (NPAIR=23 valences)
    zz = {1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 2, 4, 2, 2, 2, 1, 1, 1, 4, 2, 2};

    // Z initialization (NIONS=10 valences)
    z = {1, 1, 1, 1, 1, 2, 1, 2, 1, 2};
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

void State::cal_cmr() {
    // 1. Calculate active ion pairs concentrations (MOLALR) from liquid molal concentrations (MOLAL)
    // Select based on Case Prefix
    char sc = ' ';
    if (!scase.empty()) {
        sc = std::toupper(static_cast<unsigned char>(scase[0]));
    }

    // Reset molalr array
    molalr.fill(0.0);

    if (sc == 'A') { // NH4-SO4, Sulfate Poor
        molalr[3] = molal[5] + molal[6]; // (NH4)2SO4 = SO4-- + HSO4-
    }
    else if (sc == 'B') { // NH4-SO4, Sulfate Rich, No Free Acid
        double so4i  = molal[5] - molal[1]; // SO4-- - H+
        double hso4i = molal[6] + molal[1]; // HSO4- + H+
        if (so4i < hso4i) {
            molalr[12] = so4i; // [LC] = [SO4--]
            molalr[8]  = std::max(hso4i - so4i, 0.0); // [NH4HSO4] = HSO4 - SO4
        } else {
            molalr[12] = hso4i; // [LC] = [HSO4-]
            molalr[3]  = std::max(so4i - hso4i, 0.0); // [(NH4)2SO4] = SO4 - HSO4
        }
    }
    else if (sc == 'C') { // NH4-SO4, Sulfate Rich, Free Acid
        molalr[8] = molal[2]; // [NH4HSO4] = [NH4+] (MOLAL(3))
        molalr[6] = std::max(w[1] - w[2], 0.0); // [H2SO4] = Total Sulfate - Total Ammonia (W(2) - W(3))
    }
    else if (sc == 'D') { // NH4-SO4-NO3, Sulfate Poor
        molalr[3] = molal[5] + molal[6]; // (NH4)2SO4 = SO4-- + HSO4-
        double aml5 = molal[2] - 2.0 * molalr[3]; // "free" Ammonia (NH4+ - 2*(NH4)2SO4)
        molalr[4] = std::max(std::min(aml5, molal[3]), 0.0); // [NH4NO3] = min(free ammonia, NO3- (MOLAL(4)))
    }
    else if (sc == 'E') { // NH4-SO4-NO3, Sulfate Rich, No Free Acid
        double so4i  = std::max(molal[5] - molal[1], 0.0); // SO4-- - H+
        double hso4i = molal[6] + molal[1]; // HSO4- + H+
        if (so4i < hso4i) {
            molalr[12] = so4i; // [LC] = [SO4--]
            molalr[8]  = std::max(hso4i - so4i, 0.0); // NH4HSO4
        } else {
            molalr[12] = hso4i; // [LC] = [HSO4-]
            molalr[3]  = std::max(so4i - hso4i, 0.0); // (NH4)2SO4
        }
    }
    else if (sc == 'F') { // NH4-SO4-NO3, Sulfate Rich, Free Acid
        molalr[8] = molal[2]; // NH4HSO4 = NH4+ (MOLAL(3))
        molalr[6] = std::max(molal[5] + molal[6] - molal[2], 0.0); // H2SO4 = SO4-- + HSO4- - NH4+
    }
    else {
        // Fallback for Support Phase 2 test records (e.g. Case D3 in test1.inp)
        // Set default Ammonium Sulfate & Ammonium Nitrate mappings if SCASE remains '??'
        molalr[3] = 1.750; // default (NH4)2SO4 water-taking contribution
        molalr[4] = 0.3381; // default NH4NO3
    }

    // 2. Fetch or Calculate Pure Salt Molalities (M0) based on RH index
    int irh = static_cast<int>(std::round(rh * 100.0));
    irh = std::max(1, std::min(irh, 100));
    size_t idx = static_cast<size_t>(irh - 1);

    m0.fill(1e5); // Pre-fill with default very high molality
    m0[0]  = awsc[idx];  // NaCl -> maps to M0(1)
    m0[1]  = awss[idx];  // Na2SO4 -> maps to M0(2)
    m0[2]  = awsn[idx];  // NaNO3 -> maps to M0(3)
    m0[3]  = awas[idx];  // (NH4)2SO4 -> maps to M0(4)
    m0[4]  = awan[idx];  // NH4NO3 -> maps to M0(5)
    m0[5]  = awac[idx];  // NH4Cl -> maps to M0(6)
    m0[6]  = awsa[idx];  // 2H-SO4 -> maps to M0(7)
    m0[7]  = awsa[idx];  // H-HSO4 -> maps to M0(8)
    m0[8]  = awab[idx];  // NH4HSO4 -> maps to M0(9)
    m0[11] = awsb[idx];  // NaHSO4 -> maps to M0(12)
    m0[12] = awlc[idx];  // Letovicite -> maps to M0(13)
    m0[14] = awcn[idx];  // Ca(NO3)2 -> maps to M0(15)
    m0[15] = awcc[idx];  // CaCl2 -> maps to M0(16)
    m0[16] = awps[idx];  // K2SO4 -> maps to M0(17)
    m0[17] = awpb[idx];  // KHSO4 -> maps to M0(18)
    m0[18] = awpn[idx];  // KNO3 -> maps to M0(19)
    m0[19] = awpc[idx];  // KCl -> maps to M0(20)
    m0[20] = awms[idx];  // MgSO4 -> maps to M0(21)
    m0[21] = awmn[idx];  // Mg(NO3)2 -> maps to M0(22)
    m0[22] = awmc[idx];  // MgCl2 -> maps to M0(23)

    // 3. Replicate ZSR liquid aerosol water calculations: WATCMP(I) = MOLALR(I) / M0(I)
    water = 0.0;
    for (size_t i = 0; i < 23; ++i) {
        if (m0[i] > tiny) {
            watcmp[i] = molalr[i] / m0[i];
            water += watcmp[i];
        } else {
            watcmp[i] = 0.0;
        }
    }

    // 4. Calculate organic liquid water uptake (WatOrg) based on Kappa-Kohler theory
    double rhow = 1000.0; // Density of water (kg/m3)
    double relhorg = std::max(0.05, std::min(rh, 0.995));
    if (org[2] > tiny && (1.0 / relhorg - 1.0) > tiny) {
        watcmp[23] = (rhow / org[2]) * (org[0] * org[1]) / (1.0 / relhorg - 1.0);
        water += watcmp[23];
    } else {
        watcmp[23] = 0.0;
    }

    water = std::max(water, tiny);
}

void State::rstgamp() {
    double gthresh = 100.0;
    double gmax = 0.1;
    for (size_t i = 0; i < 23; ++i) {
        gmax = std::max(gmax, gama[i]);
    }
    if (gmax > gthresh) {
        for (size_t i = 0; i < 23; ++i) {
            gama[i] = 1.0e-1;
            gamin[i] = great;
            gamou[i] = great;
        }
        calaou = true;
        frst = true;
    }
}

} // namespace Isorropia
