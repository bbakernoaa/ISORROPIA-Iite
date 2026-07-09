#ifndef ISORROPIA_C_API_H
#define ISORROPIA_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Plain C-layout input structure for passing simulation parameters.
 * 
 * Sized and ordered to align perfectly across compilers with standard C layouts.
 */
typedef struct {
    double w[8];      ///< Total concentrations (µmol/m3 of air or ug/m3).
    double org[3];    ///< Organic parameters (concentration, k_org, density).
    double waer[8];   ///< Aerosol-phase concentrations.
    double temp;      ///< Ambient temperature (Kelvin).
    double rh;        ///< Relative Humidity (fraction [0.0 - 1.0]).
    int iprob;        ///< Problem formulation (0 = Forward, 1 = Reverse).
    int nadj;         ///< Mass adjustment toggle.
} IsorropiaInput;

/**
 * @brief Plain C-layout structure representing the complete thermodynamic state.
 * 
 * All elements map directly to members of the thread-safe C++ Isorropia::State structure.
 * This can be safely cast or read by C/C++ or Fortran host model programs.
 */
typedef struct {
    double temp;
    double rh;
    double w[8];
    double waer[8];
    double org[3];

    double molal[10];
    double molalr[23];
    double gama[23];
    double zz[23];
    double z[10];
    double gamou[23];
    double gamin[23];
    double m0[23];
    double gasaq[3];
    int actmod;
    double epsact;
    double coh;
    double chno3;
    double chcl;
    double water;
    double ionic;
    double watcmp[24];
    int frst;
    int calain;
    int calaou;
    int dryf;

    double ch2so4, cnh42s4, cnh4hs4, cnacl, cna2so4, cnano3, cnh4no3, cnh4cl, cnahso4, clc;
    double ccaso4, ccano32, ccacl2, ck2so4, ckhso4, ckno3, ckcl, cmgso4, cmgno32, cmgcl2;
    double gnh3, ghno3, ghcl;

    int num_errors;
} IsorropiaState;

/**
 * @brief Top-level execution entry point for C, C++, and Fortran binders.
 * 
 * Executes the atmospheric aerosol thermodynamics solver using thread-safe, 
 * standard C linkages (extern "C").
 * 
 * @param input Pointer to standard C plain layout Input parameters.
 * @param state Pointer to standard C plain layout output State structure.
 */
void isorropia_solve_c(const IsorropiaInput* input, IsorropiaState* state);

#ifdef __cplusplus
}
#endif

#endif // ISORROPIA_C_API_H
