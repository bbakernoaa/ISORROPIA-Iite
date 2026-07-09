#ifndef ISORROPIA_SOLVER_HPP
#define ISORROPIA_SOLVER_HPP

#include <array>
#include <string>
#include <string_view>

namespace Isorropia {

/**
 * @brief Species components indices mapped to 0-based array elements.
 * 
 * This enum maps components to their index in the 0-indexed arrays, replacing
 * the legacy 1-indexed Fortran array layout for components (W).
 */
enum class Component : size_t {
    Na    = 0,  ///< Sodium (Na+) ion / total sodium
    H2SO4 = 1,  ///< Sulfuric acid (H2SO4) / total sulfate
    NH3   = 2,  ///< Ammonia (NH3) / total ammonia
    HNO3  = 3,  ///< Nitric acid (HNO3) / total nitrate
    HCl   = 4,  ///< Hydrochloric acid (HCl) / total chloride
    Ca    = 5,  ///< Calcium (Ca2+) ion / total calcium
    K     = 6,  ///< Potassium (K+) ion / total potassium
    Mg    = 7   ///< Magnesium (Mg2+) ion / total magnesium
};

/**
 * @brief Liquid aerosol ion species indices mapped to 0-based elements.
 * 
 * Maps to Fortran 'MOLAL' array in COMMON /IONS/ (1-indexed MOLAL(NIONS)).
 */
enum class Ions : size_t {
    Na   = 0,  ///< Liquid Sodium ion [Na+]
    H    = 1,  ///< Liquid Hydrogen ion [H+]
    NH4  = 2,  ///< Liquid Ammonium ion [NH4+]
    NO3  = 3,  ///< Liquid Nitrate ion [NO3-]
    Cl   = 4,  ///< Liquid Chloride ion [Cl-]
    SO4  = 5,  ///< Liquid Sulfate ion [SO4^2-]
    HSO4 = 6,  ///< Liquid Bisulfate ion [HSO4-]
    Ca   = 7,  ///< Liquid Calcium ion [Ca^2+]
    K    = 8,  ///< Liquid Potassium ion [K+]
    Mg   = 9   ///< Liquid Magnesium ion [Mg^2+]
};

/**
 * @brief Input parameters for a single simulation cell.
 * 
 * Encapsulates all variables originally shared via COMMON /INPT/ in 'isrpia.inc'.
 */
struct Input {
    /**
     * @brief Total concentrations of species in the system.
     * 
     * Unit: µmol/m3 of air (when input unit is 0) or ug/m3 of air (when input unit is 1).
     * Maps to legacy Fortran 'W(NCOMP)' array (W(8)).
     */
    std::array<double, 8> w = {0.0};

    /**
     * @brief Organic species concentrations and properties.
     * 
     * Unit: ug/m3 of air (for concentration), fraction (for organic-water uptake k_org), and kg/m3 (for density).
     * Maps to legacy Fortran 'ORG(NORG)' array (ORG(3)).
     * - org[0]: Concentration of organics (ug/m3)
     * - org[1]: Organic hygroscopicity parameter (k_org)
     * - org[2]: Density of organic phase (rho_org, kg/m3)
     */
    std::array<double, 3> org = {0.0};

    /**
     * @brief Concentration of aerosol-phase components.
     * 
     * Unit: µmol/m3 of air or ug/m3.
     * Maps to legacy Fortran 'WAER(NCOMP)' array.
     */
    std::array<double, 8> waer = {0.0};

    /**
     * @brief Ambient temperature.
     * 
     * Unit: Kelvin (K).
     * Maps to legacy Fortran variable 'TEMP'.
     */
    double temp = 298.15;

    /**
     * @brief Ambient Relative Humidity.
     * 
     * Unit: Fraction [0.0 to 1.0] (e.g. 0.80 for 80% RH).
     * Maps to legacy Fortran variable 'RH'.
     */
    double rh = 0.0;

    /**
     * @brief Problem type / formulation selector.
     * 
     * - 0: Forward problem (Total concentrations are specified as inputs).
     * - 1: Reverse problem (Only aerosol-phase concentrations are specified as inputs).
     * Maps to legacy Fortran variable 'IPROB'.
     */
    int iprob = 0;

    /**
     * @brief Mass balance adjustment switch.
     * 
     * - 0: No adjustment.
     * - 1: Adjust concentrations for mass balance.
     * Maps to legacy Fortran variable 'NADJ'.
     */
    int nadj = 0;
};

/**
 * @brief Single entry in the thread-local non-throwing error stack.
 */
struct ErrorEntry {
    int code = 0;              ///< Numeric error code returned by the solver.
    std::string message;       ///< Descriptive error message explaining the issue.
};

/**
 * @brief Comprehensive internal state of a simulation cell.
 * 
 * Fully encapsulates all legacy Fortran COMMON blocks to ensure total thread safety.
 * This class contains no global data, meaning multiple threads can execute solvers simultaneously on separate states.
 */
struct State {
    //=======================================================================
    // COMMON /INPT/ equivalents (Meteorological state and inputs)
    //=======================================================================
    double temp = 298.15;  ///< Ambient temperature (Kelvin). Maps to Fortran 'TEMP' in COMMON /INPT/.
    double rh = 0.0;       ///< Ambient relative humidity (fraction [0.0 - 1.0]). Maps to Fortran 'RH' in COMMON /INPT/.
    std::array<double, 8> w = {0.0};     ///< Total concentrations copy. Maps to Fortran 'W(NCOMP)' in COMMON /INPT/.
    std::array<double, 8> waer = {0.0};  ///< Aerosol concentrations copy. Maps to Fortran 'WAER(NCOMP)' in COMMON /INPT/.
    std::array<double, 3> org = {0.0};   ///< Organic concentrations copy. Maps to Fortran 'ORG(NORG)' in COMMON /INPT/.

    //=======================================================================
    // COMMON /IONS/ equivalents (Liquid aerosol phase properties)
    //=======================================================================
    
    /**
     * @brief Molal concentrations of liquid ions.
     * 
     * Unit: mol/kg of water.
     * Maps to legacy Fortran 'MOLAL(NIONS)' (MOLAL(10)).
     */
    std::array<double, 10> molal = {0.0};

    /**
     * @brief Molal concentrations of active ion pairs.
     * 
     * Maps to legacy Fortran 'MOLALR(NPAIR)' (MOLALR(23)).
     */
    std::array<double, 23> molalr = {0.0};

    /**
     * @brief Mean molal activity coefficients of active ion pairs.
     * 
     * Maps to legacy Fortran 'GAMA(NPAIR)' (GAMA(23)).
     */
    std::array<double, 23> gama = {0.0};

    /**
     * @brief Charges of active ion pairs.
     * 
     * Maps to legacy Fortran 'ZZ(NPAIR)' (ZZ(23)).
     */
    std::array<double, 23> zz = {0.0};

    /**
     * @brief Charges of individual ion species.
     * 
     * Maps to legacy Fortran 'Z(NIONS)' (Z(10)).
     */
    std::array<double, 10> z = {0.0};

    /**
     * @brief Activity coefficients in the outer iteration shell.
     * 
     * Maps to legacy Fortran 'GAMOU(NPAIR)'.
     */
    std::array<double, 23> gamou = {0.0};

    /**
     * @brief Activity coefficients in the inner iteration shell.
     * 
     * Maps to legacy Fortran 'GAMIN(NPAIR)'.
     */
    std::array<double, 23> gamin = {0.0};

    /**
     * @brief Reference molality constants for active ion pairs.
     * 
     * Maps to legacy Fortran 'M0(NPAIR)'.
     */
    std::array<double, 23> m0 = {0.0};

    /**
     * @brief Dissolved gaseous species in liquid aerosol phase.
     * 
     * Unit: mol/m3 of air.
     * Maps to legacy Fortran 'GASAQ(NGASAQ)' (GASAQ(3)).
     * - gasaq[0]: Aquated Ammonia (NH3_aq)
     * - gasaq[1]: Aquated Nitric Acid (HNO3_aq)
     * - gasaq[2]: Aquated Hydrochloric Acid (HCl_aq)
     */
    std::array<double, 3> gasaq = {0.0};

    int actmod = 0;       ///< Activity model coefficient selector. Maps to Fortran 'ACTMOD'.
    double epsact = 0.0;  ///< Activity calculation precision tolerance. Maps to Fortran 'EPSACT'.
    double coh = 0.0;     ///< Liquid phase Hydroxide (OH-) concentration. Maps to Fortran 'COH'.
    double chno3 = 0.0;   ///< Dissolved liquid phase HNO3. Maps to Fortran 'CHNO3'.
    double chcl = 0.0;    ///< Dissolved liquid phase HCl. Maps to Fortran 'CHCL'.
    
    /**
     * @brief Total liquid water content of the aerosol phase.
     * 
     * Unit: kg/m3 of air.
     * Maps to legacy Fortran variable 'WATER'.
     */
    double water = 0.0;
    
    double ionic = 0.0;   ///< Total ionic strength of liquid aerosol. Maps to Fortran 'IONIC'.

    /**
     * @brief Individual component contributions of salts/organics to liquid water content.
     * 
     * Unit: kg/m3 of air.
     * Maps to legacy Fortran 'WATCMP(NPAIR+1)' (WATCMP(24)).
     */
    std::array<double, 24> watcmp = {0.0};

    bool frst = true;      ///< First execution run initialization flag. Maps to Fortran 'FRST'.
    bool calain = false;   ///< Inner activity coefficient calculation flag. Maps to Fortran 'CALAIN'.
    bool calaou = false;   ///< Outer activity coefficient calculation flag. Maps to Fortran 'CALAOU'.
    bool dryf = false;     ///< True if aerosol phase is completely dry. Maps to Fortran 'DRYF'.

    //=======================================================================
    // COMMON /SALT/ equivalents (Solid aerosol phase species concentrations)
    // Unit: µmol/m3 or mol/m3 (depending on solver output mapping)
    //=======================================================================
    double ch2so4 = 0.0;   ///< Solid Sulfuric acid (H2SO4). Maps to Fortran 'CH2SO4'.
    double cnh42s4 = 0.0;  ///< Solid Ammonium Sulfate ((NH4)2SO4). Maps to Fortran 'CNH42S4'.
    double cnh4hs4 = 0.0;  ///< Solid Ammonium Bisulfate (NH4HSO4). Maps to Fortran 'CNH4HS4'.
    double cnacl = 0.0;    ///< Solid Sodium Chloride (NaCl). Maps to Fortran 'CNACL'.
    double cna2so4 = 0.0;  ///< Solid Sodium Sulfate (Na2SO4). Maps to Fortran 'CNA2SO4'.
    double cnano3 = 0.0;   ///< Solid Sodium Nitrate (NaNO3). Maps to Fortran 'CNANO3'.
    double cnh4no3 = 0.0;  ///< Solid Ammonium Nitrate (NH4NO3). Maps to Fortran 'CNH4NO3'.
    double cnh4cl = 0.0;   ///< Solid Ammonium Chloride (NH4Cl). Maps to Fortran 'CNH4CL'.
    double cnahso4 = 0.0;  ///< Solid Sodium Bisulfate (NaHSO4). Maps to Fortran 'CNAHSO4'.
    double clc = 0.0;      ///< Solid Letovicite ((NH4)3H(SO4)2). Maps to Fortran 'CLC'.
    double ccaso4 = 0.0;   ///< Solid Calcium Sulfate (CaSO4). Maps to Fortran 'CCASO4'.
    double ccano32 = 0.0;  ///< Solid Calcium Nitrate (Ca(NO3)2). Maps to Fortran 'CCANO32'.
    double ccacl2 = 0.0;   ///< Solid Calcium Chloride (CaCl2). Maps to Fortran 'CCACL2'.
    double ck2so4 = 0.0;   ///< Solid Potassium Sulfate (K2SO4). Maps to Fortran 'CK2SO4'.
    double ckhso4 = 0.0;   ///< Solid Potassium Bisulfate (KHSO4). Maps to Fortran 'CKHSO4'.
    double ckno3 = 0.0;    ///< Solid Potassium Nitrate (KNO3). Maps to Fortran 'CKNO3'.
    double ckcl = 0.0;     ///< Solid Potassium Chloride (KCl). Maps to Fortran 'CKCL'.
    double cmgso4 = 0.0;   ///< Solid Magnesium Sulfate (MgSO4). Maps to Fortran 'CMGSO4'.
    double cmgno32 = 0.0;  ///< Solid Magnesium Nitrate (Mg(NO3)2). Maps to Fortran 'CMGNO32'.
    double cmgcl2 = 0.0;   ///< Solid Magnesium Chloride (MgCl2). Maps to Fortran 'CMGCL2'.

    //=======================================================================
    // COMMON /GAS/ equivalents (Gas phase species concentrations)
    // Unit: mol/m3 of air.
    //=======================================================================
    double gnh3 = 0.0;     ///< Gaseous Ammonia (NH3). Maps to Fortran 'GNH3'.
    double ghno3 = 0.0;    ///< Gaseous Nitric Acid (HNO3). Maps to Fortran 'GHNO3'.
    double ghcl = 0.0;     ///< Gaseous Hydrochloric Acid (HCl). Maps to Fortran 'GHCL'.

    //=======================================================================
    // COMMON /DRH / equivalents (Unicomponent Deliquescence Relative Humidities)
    //=======================================================================
    double drh2so4 = 0.0;   ///< DRH of Sulfuric acid (H2SO4). Maps to Fortran 'DRH2SO4'.
    double drnh42s4 = 0.0;  ///< DRH of Ammonium Sulfate ((NH4)2SO4). Maps to Fortran 'DRNH42S4'.
    double drnahso4 = 0.0;  ///< DRH of Sodium Bisulfate (NaHSO4). Maps to Fortran 'DRNAHSO4'.
    double drnacl = 0.0;    ///< DRH of Sodium Chloride (NaCl). Maps to Fortran 'DRNACL'.
    double drnano3 = 0.0;   ///< DRH of Sodium Nitrate (NaNO3). Maps to Fortran 'DRNANO3'.
    double drna2so4 = 0.0;  ///< DRH of Sodium Sulfate (Na2SO4). Maps to Fortran 'DRNA2SO4'.
    double drnh4hs4 = 0.0;  ///< DRH of Ammonium Bisulfate (NH4HSO4). Maps to Fortran 'DRNH4HS4'.
    double drlc = 0.0;      ///< DRH of Letovicite ((NH4)3H(SO4)2). Maps to Fortran 'DRLC'.
    double drnh4no3 = 0.0;  ///< DRH of Ammonium Nitrate (NH4NO3). Maps to Fortran 'DRNH4NO3'.
    double drnh4cl = 0.0;   ///< DRH of Ammonium Chloride (NH4Cl). Maps to Fortran 'DRNH4CL'.
    double drcaso4 = 0.0;   ///< DRH of Calcium Sulfate (CaSO4). Maps to Fortran 'DRCASO4'.
    double drcano32 = 0.0;  ///< DRH of Calcium Nitrate (Ca(NO3)2). Maps to Fortran 'DRCANO32'.
    double drcacl2 = 0.0;   ///< DRH of Calcium Chloride (CaCl2). Maps to Fortran 'DRCACL2'.
    double drk2so4 = 0.0;   ///< DRH of Potassium Sulfate (K2SO4). Maps to Fortran 'DRK2SO4'.
    double drkhso4 = 0.0;   ///< DRH of Potassium Bisulfate (KHSO4). Maps to Fortran 'DRKHSO4'.
    double drkno3 = 0.0;    ///< DRH of Potassium Nitrate (KNO3). Maps to Fortran 'DRKNO3'.
    double drkcl = 0.0;     ///< DRH of Potassium Chloride (KCl). Maps to Fortran 'DRKCL'.
    double drmgso4 = 0.0;   ///< DRH of Magnesium Sulfate (MgSO4). Maps to Fortran 'DRMGSO4'.
    double drmgno32 = 0.0;  ///< DRH of Magnesium Nitrate (Mg(NO3)2). Maps to Fortran 'DRMGNO32'.
    double drmgcl2 = 0.0;   ///< DRH of Magnesium Chloride (MgCl2). Maps to Fortran 'DRMGCL2'.

    //=======================================================================
    // COMMON /MDRH/ & /MDRH2/ equivalents (Mutual Deliquescence Relative Humidities)
    //=======================================================================
    double drmlcab = 0.0, drmlcas = 0.0, drmasan = 0.0, drmg1 = 0.0, drmg2 = 0.0;
    double drmg3 = 0.0, drmh1 = 0.0, drmh2 = 0.0, drmi1 = 0.0, drmi2 = 0.0;
    double drmi3 = 0.0, drmq1 = 0.0, drmr1 = 0.0, drmr2 = 0.0, drmr3 = 0.0;
    double drmr4 = 0.0, drmr5 = 0.0, drmr6 = 0.0, drmr7 = 0.0, drmr8 = 0.0;
    double drmr9 = 0.0, drmr10 = 0.0, drmr11 = 0.0, drmr12 = 0.0, drmr13 = 0.0;
    int wftyp = 0;          ///< Water formulation type. Maps to Fortran 'WFTYP'.

    double drmo1 = 0.0, drmo2 = 0.0, drmo3 = 0.0, drml1 = 0.0, drml2 = 0.0;
    double drml3 = 0.0, drmmm1 = 0.0, drmmm2 = 0.0, drmp1 = 0.0, drmp2 = 0.0;
    double drmp3 = 0.0, drmp4 = 0.0, drmp5 = 0.0, drmv1 = 0.0;

    //=======================================================================
    // COMMON /EQUK/ equivalents (Chemical Equilibrium Constants)
    //=======================================================================
    double xk1  = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK1'.
    double xk2  = 0.0;  ///< Equilibrium constant (derived: xk21*xk22). Maps to Fortran 'XK2'.
    double xk3  = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK3'.
    double xk4  = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK4'.
    double xk5  = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK5'.
    double xk6  = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK6'.
    double xk7  = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK7'.
    double xk8  = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK8'.
    double xk9  = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK9'.
    double xk10 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK10'.
    double xk11 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK11'.
    double xk12 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK12'.
    double xk13 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK13'.
    double xk14 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK14'.
    double xkw  = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XKW'.
    double xk21 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK21'.
    double xk22 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK22'.
    double xk31 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK31'.
    double xk32 = 0.0;  ///< Equilibrium constant (derived: xk3/xk31). Maps to Fortran 'XK32'.
    double xk41 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK41'.
    double xk42 = 0.0;  ///< Equilibrium constant (derived: xk4/xk41). Maps to Fortran 'XK42'.
    double xk15 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK15'.
    double xk16 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK16'.
    double xk17 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK17'.
    double xk18 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK18'.
    double xk19 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK19'.
    double xk20 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK20'.
    double xk23 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK23'.
    double xk24 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK24'.
    double xk25 = 0.0;  ///< Equilibrium constant. Maps to Fortran 'XK25'.

    //=======================================================================
    // COMMON /CASE/ equivalents (Simulation Case and Regimes)
    //=======================================================================
    std::string scase = "??";  ///< Solution regime/case description. Maps to Fortran 'SCASE' in COMMON /CASE/.

    //=======================================================================
    // COMMON /ZSR/ equivalents (Water activities arrays)
    //=======================================================================
    std::array<double, 100> awas = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWAS(NZSR)'.
    std::array<double, 100> awss = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWSS(NZSR)'.
    std::array<double, 100> awac = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWAC(NZSR)'.
    std::array<double, 100> awsc = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWSC(NZSR)'.
    std::array<double, 100> awan = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWAN(NZSR)'.
    std::array<double, 100> awsn = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWSN(NZSR)'.
    std::array<double, 100> awsb = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWSB(NZSR)'.
    std::array<double, 100> awab = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWAB(NZSR)'.
    std::array<double, 100> awsa = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWSA(NZSR)'.
    std::array<double, 100> awlc = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWLC(NZSR)'.
    std::array<double, 100> awcs = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWCS(NZSR)'.
    std::array<double, 100> awcn = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWCN(NZSR)'.
    std::array<double, 100> awcc = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWCC(NZSR)'.
    std::array<double, 100> awps = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWPS(NZSR)'.
    std::array<double, 100> awpb = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWPB(NZSR)'.
    std::array<double, 100> awpn = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWPN(NZSR)'.
    std::array<double, 100> awpc = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWPC(NZSR)'.
    std::array<double, 100> awms = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWMS(NZSR)'.
    std::array<double, 100> awmn = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWMN(NZSR)'.
    std::array<double, 100> awmc = {0.0};  ///< ZSR Water activity grid. Maps to Fortran 'AWMC(NZSR)'.

    //=======================================================================
    // Thread-Safe Non-Throwing Error Stack
    //=======================================================================
    bool stack_overflow = false;                ///< Set to true if more than 25 errors are pushed.
    size_t num_errors = 0;                      ///< Active number of errors currently on the stack.
    std::array<ErrorEntry, 25> error_stack;     ///< Pre-allocated stack to avoid dynamic allocation overhead.

    /**
     * @brief Pushes a diagnostic error onto the local stack without throwing an exception.
     * 
     * Keeps execution thread-safe and extremely fast.
     * @param code The numeric error identifier.
     * @param message Explanatory text for the error.
     */
    void push_error(int code, std::string_view message) {
        if (num_errors >= error_stack.size()) {
            stack_overflow = true;
            return;
        }
        error_stack[num_errors] = {code, std::string(message)};
        num_errors++;
    }

    /**
     * @brief Clears all errors from the thread-local stack.
     */
    void clear_errors() {
        num_errors = 0;
        stack_overflow = false;
    }

    //=======================================================================
    // Physical Constants and Molecular Weights (BLOCK DATA BLKISO equivalents)
    //=======================================================================
    double r = 82.0567e-6;  ///< Gas constant. Maps to Fortran 'R' (m3 atm / mol K).
    double tiny = 1e-20;    ///< Small threshold. Maps to Fortran 'TINY'.
    double tiny2 = 1e-11;   ///< Small threshold 2. Maps to Fortran 'TINY2'.
    double great = 1e10;    ///< Large threshold. Maps to Fortran 'GREAT'.
    double zero = 0.0;      ///< Constant zero. Maps to Fortran 'ZERO'.
    double one = 1.0;       ///< Constant one. Maps to Fortran 'ONE'.

    // Molecular weights
    std::array<double, 10> imw = {0.0};  ///< Molecular weights of 10 ions (g/mol). Maps to Fortran 'IMW(NIONS)'.
                                         ///< 0: Na+ (23.0), 1: H+ (1.0), 2: NH4+ (18.0), 3: NO3- (62.0), 4: Cl- (35.5),
                                         ///< 5: SO4^2- (96.0), 6: HSO4- (97.0), 7: Ca^2+ (40.1), 8: K+ (39.1), 9: Mg^2+ (24.3)
    
    std::array<double, 8> wmw = {0.0};   ///< Molecular weights of 8 components (g/mol). Maps to Fortran 'WMW(NCOMP)'.
                                         ///< 0: Na (23.0), 1: H2SO4 (98.0), 2: NH3 (17.0), 3: HNO3 (63.0), 4: HCl (36.5),
                                         ///< 5: Ca (40.1), 6: K (39.1), 7: Mg (24.3)

    std::array<double, 23> smw = {0.0};  ///< Molecular weights of 23 active salt pairs (g/mol). Maps to Fortran 'SMW(NPAIR)'.
                                         ///< 0: NaNO3 (85.0), 1: NH4NO3 (80.0), 2: NaCl (58.5), 3: NH4Cl (53.5),
                                         ///< 4: Na2SO4 (142.0), 5: (NH4)2SO4 (132.0), 6: NaHSO4 (120.0), 7: NH4HSO4 (115.0),
                                         ///< 8: H2SO4_aq (120.0), 9: Ca(NO3)2 (136.0), 10: CaSO4 (164.0), 11: CaCl2 (111.0),
                                         ///< 12: K2SO4 (174.0), 13: KHSO4 (136.0), 14: KNO3 (101.0), 15: KCl (74.5),
                                         ///< 16: MgSO4 (120.0), 17: Mg(NO3)2 (148.0), 18: MgCl2 (95.0), 19: H2O (18.0),
                                         ///< 20-22: placeholders (0.0)

    /**
     * @brief Initializes molecular weights and fundamental physical constants.
     * 
     * Replicates block data initialization in Fortran 'BLOCK DATA BLKISO'.
     */
    void initialize_constants();

    /**
     * @brief Initializes ZSR water activity lookup tables.
     * 
     * Replicates pure salt data grids from Fortran 'BLOCK DATA BLKISO'.
     */
    void initialize_water_activities();

    /**
     * @brief Initializes unicomponent and mutual Deliquescence Relative Humidities.
     * 
     * Replicates DRH calculations and temperature dependency formulas from Fortran INIT subroutines.
     */
    void initialize_drh();

    /**
     * @brief Calculates temperature-dependent equilibrium constants.
     * 
     * Replicates calculations and corrections from Fortran subroutines (van 't Hoff equation).
     */
    void calculate_equilibrium_constants();

    /**
     * @brief Computes binary activity coefficients using pre-tabulated Kusik-Meissner grids.
     * 
     * Replaces Fortran 'SUBROUTINE KMTAB'.
     * @param ionic_strength Total ionic strength of the solution.
     * @param temp_k Current temperature of the cell in Kelvin (K).
     * @param g0 Output array populated with 23 binary activity coefficients.
     */
    void km_tab(double ionic_strength, double temp_k, std::array<double, 23>& g0);

    /**
     * @brief Computes multicomponent activity coefficients for pure Case 1 systems.
     * 
     * Replaces Fortran 'SUBROUTINE CALCACT1'.
     */
    void cal_act1();

    /**
     * @brief Computes dynamic multi-species liquid water content of the aerosol using the ZSR relation.
     * 
     * Replaces Fortran 'SUBROUTINE CALCMR'.
     */
    void cal_cmr();
};

/**
 * @brief Thread-safe implementation of the ISORROPIA-Lite aerosol thermodynamics solver.
 */
class Solver {
public:
    Solver();

    /**
     * @brief Solves the thermodynamic equilibrium for the given inputs and state.
     * 
     * Replaces the legacy Fortran 'SUBROUTINE ISOROPIA'.
     * @param input The chemical and meteorological inputs for the grid cell.
     * @param state The local state structure where results are calculated and stored.
     */
    void solve(const Input& input, State& state);

private:
    /**
     * @brief Forward solver for NH4-SO4-H2O systems (Case 1).
     * 
     * Maps to legacy Fortran 'SUBROUTINE ISRP1F' in 'isofwd.f'.
     */
    void isrp1f(const Input& input, State& state);

    /**
     * @brief Forward solver for NH4-SO4-NO3-H2O systems (Case 2).
     * 
     * Maps to legacy Fortran 'SUBROUTINE ISRP2F' in 'isofwd.f'.
     */
    void isrp2f(const Input& input, State& state);

    /**
     * @brief Forward solver for Na-NH4-SO4-NO3-Cl-H2O systems (Case 3).
     * 
     * Maps to legacy Fortran 'SUBROUTINE ISRP3F' in 'isofwd.f'.
     */
    void isrp3f(const Input& input, State& state);

    /**
     * @brief Forward solver for Na-NH4-SO4-NO3-Cl-Ca-K-Mg-H2O crustal systems (Case 4).
     * 
     * Maps to legacy Fortran 'SUBROUTINE ISRP4F' in 'isofwd.f'.
     */
    void isrp4f(const Input& input, State& state);
};

} // namespace Isorropia

#endif // ISORROPIA_SOLVER_HPP
