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
    
    float ionic = 0.0f;   ///< Total ionic strength of liquid aerosol. Maps to Fortran 'IONIC'.

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
};

} // namespace Isorropia

#endif // ISORROPIA_SOLVER_HPP
