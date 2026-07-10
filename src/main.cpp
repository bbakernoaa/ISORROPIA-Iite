#include "Isorropia/Solver.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>

struct InputRecord {
    double na = 0.0;
    double so4 = 0.0;
    double nh3 = 0.0;
    double no3 = 0.0;
    double cl = 0.0;
    double ca = 0.0;
    double k = 0.0;
    double mg = 0.0;
    double org = 0.0;
    double k_org = 0.0;
    double rho_org = 0.0;
    double rh = 0.0;
    double temp = 0.0;
};

// Simple helper to trim strings
std::string trim_str(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n\"");
    if (std::string::npos == first) {
        return "";
    }
    size_t last = str.find_last_not_of(" \t\r\n\"");
    return str.substr(first, (last - first + 1));
}

int main(int argc, char* argv[]) {
    std::string input_file = "test1.inp";
    if (argc > 1) {
        input_file = argv[1];
    }

    std::cout << "=== ISORROPIA-Lite C++ CLI Driver ===" << std::endl;
    std::cout << "Reading input file: " << input_file << std::endl;

    std::ifstream infile(input_file);
    if (!infile.is_open()) {
        // Try looking in sibling dirs
        infile.open("isolite1_0_src/" + input_file);
        if (!infile.is_open()) {
            infile.open("../isolite1_0_src/" + input_file);
            if (!infile.is_open()) {
                infile.open("ISORROPIALite_Executable_Manual_Papers/" + input_file);
                if (!infile.is_open()) {
                    infile.open("../ISORROPIALite_Executable_Manual_Papers/" + input_file);
                }
            }
        }
    }

    if (!infile.is_open()) {
        std::cerr << "Error: Could not open input file " << input_file << std::endl;
        return 1;
    }

    std::string line;
    bool headers_done = false;
    std::vector<InputRecord> records;

    // Control values parsed from the header of the .inp file, matching the
    // legacy Fortran main.f reading order:
    //   1st numeric line  -> input units (INUNIT): 0 = umol/m3, 1 = ug/m3
    //   2nd numeric line  -> problem type (IPROB): 0 = forward, 1 = reverse
    //                        (the trailing phase-state field is intentionally
    //                         ignored; ISORROPIA-Lite is metastable-only, just
    //                         like the reference Fortran binary).
    std::vector<int> control_values;

    while (std::getline(infile, line)) {
        std::string trimmed = trim_str(line);
        if (trimmed.empty() || trimmed[0] == 'C' || trimmed[0] == '#' || trimmed[0] == '*') continue;

        if (!headers_done) {
            // Locate species column header row
            if (trimmed.find("Na") != std::string::npos && trimmed.find("SO4") != std::string::npos) {
                headers_done = true;
                continue;
            }
            // Otherwise, capture leading numeric control tokens (units, iprob).
            std::stringstream cs(trimmed);
            double ctrl_val;
            if (cs >> ctrl_val) {
                control_values.push_back(static_cast<int>(std::lround(ctrl_val)));
            }
            continue;
        }

        std::stringstream ss(trimmed);
        InputRecord rec;
        if (ss >> rec.na >> rec.so4 >> rec.nh3 >> rec.no3 >> rec.cl >> rec.ca >> rec.k >> rec.mg >> rec.org >> rec.k_org >> rec.rho_org >> rec.rh >> rec.temp) {
            records.push_back(rec);
        }
    }

    // Resolve control values (defaults match legacy behavior: ug/m3, forward).
    int input_units = control_values.size() > 0 ? control_values[0] : 1;
    int iprob       = control_values.size() > 1 ? control_values[1] : 0;

    std::cout << "Parsed " << records.size() << " records from input file." << std::endl;
    std::cout << "Input units: " << (input_units == 0 ? "umol/m3" : "ug/m3")
              << " | Problem type: " << (iprob == 1 ? "REVERSE" : "FORWARD")
              << " | Aerosol state: METASTABLE (liquid only)" << std::endl;

    // Output filename
    std::string base_name = input_file;
    size_t last_slash = base_name.find_last_of("/\\");
    if (last_slash != std::string::npos) {
        base_name = base_name.substr(last_slash + 1);
    }
    size_t dot_pos = base_name.find_last_of('.');
    if (dot_pos != std::string::npos) {
        base_name = base_name.substr(0, dot_pos);
    }
    std::string out_file = base_name + "_cpp.txt";

    std::ofstream outfile(out_file);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file " << out_file << std::endl;
        return 1;
    }

    Isorropia::Solver solver;
    const std::vector<double> wmw = {23.0, 98.0, 17.0, 63.0, 36.5, 40.1, 39.1, 24.3}; // Component Molecular Weights

    for (size_t r = 0; r < records.size(); ++r) {
        const auto& rec = records[r];

        Isorropia::Input input;
        // Convert total/aerosol inputs to mol/m3, honoring the input-unit flag
        // (mirrors Fortran INPDAT: umol/m3 -> *1e-6; ug/m3 -> /MW *1e-6).
        const double raw[8] = {rec.na, rec.so4, rec.nh3, rec.no3, rec.cl, rec.ca, rec.k, rec.mg};
        for (size_t i = 0; i < 8; ++i) {
            double mol_per_m3 = (input_units == 0)
                ? raw[i] * 1e-6          // umol/m3 -> mol/m3
                : (raw[i] / wmw[i]) * 1e-6; // ug/m3   -> mol/m3
            input.w[i] = std::max(mol_per_m3, 0.0);
        }

        // Convert Organic concentrations to kilograms/m3 (ug/m3 to kg/m3) and density to kg/m3
        input.org[0] = rec.org * 1e-9;
        input.org[1] = rec.k_org;
        input.org[2] = rec.rho_org * 1e3;

        input.rh = rec.rh;
        input.temp = rec.temp;

        // Problem type comes from the .inp control line (not the filename).
        if (iprob == 1) {
            input.iprob = 1; // Reverse problem: inputs are aerosol-phase totals
            std::copy(input.w.begin(), input.w.end(), input.waer.begin());
            std::fill(input.w.begin(), input.w.end(), 0.0);
        } else {
            input.iprob = 0; // Forward problem
        }
        input.nadj = 1;


        Isorropia::State state;
        solver.solve(input, state);

        // Dynamically compute correct outputs in ug/m3 to match test1.txt reference columns
        double water_val = state.water * 1e9; // water in ug/m3 (since state.water is in kg/m3)

        double h_print     = state.molal[1] * 1.0 * 1e6;  // mol/m3 * MW * 1e6 -> ug/m3
        double nh4_print   = state.molal[2] * 18.0 * 1e6; // NH4+ (MW=18)
        double no3_print   = state.molal[3] * 62.0 * 1e6; // NO3- (MW=62)
        double so4_print   = state.molal[5] * 96.0 * 1e6; // SO4-- (MW=96)
        double hso4_print  = state.molal[6] * 97.0 * 1e6; // HSO4- (MW=97)
        double na_print    = state.molal[0] * 23.0 * 1e6; // Na+ (MW=23)
        double cl_print    = state.molal[4] * 35.5 * 1e6; // Cl- (MW=35.5)
        double ca_print    = state.molal[7] * 40.1 * 1e6; // Ca++ (MW=40.1)
        double k_print     = state.molal[8] * 39.1 * 1e6; // K+ (MW=39.1)
        double mg_print    = state.molal[9] * 24.3 * 1e6; // Mg++ (MW=24.3)

        double gnh3_print  = state.gnh3 * 17.0 * 1e6; // Gaseous Ammonia (MW=17) -> ug/m3
        double ghno3_print = state.ghno3 * 63.0 * 1e6; // Gaseous Nitric Acid (MW=63) -> ug/m3
        double ghcl_print  = state.ghcl * 36.5 * 1e6; // Gaseous HCl (MW=36.5) -> ug/m3

        double ph_print    = (state.water > 1e-20 && state.molal[1] > 1e-30) ? -std::log10(state.molal[1] / state.water) : 7.0;
        double ionic_print = state.ionic;

        double wat_nh42so4 = state.watcmp[3] * 1e9;
        double wat_nh4no3  = state.watcmp[4] * 1e9;
        double wat_org     = state.watcmp[23] * 1e9;

        outfile << "Record " << r + 1 << " Solution:" << std::endl;
        outfile << " SCASE " << state.scase << "\n";
        outfile << " [WATER ] " << std::scientific << std::setprecision(3) << water_val << "\n";
        outfile << " [H+ ] " << std::scientific << std::setprecision(3) << h_print << "\n";
        outfile << " [Na+ ] " << std::scientific << std::setprecision(3) << na_print << "\n";
        outfile << " [NH4+ ] " << std::scientific << std::setprecision(3) << nh4_print << "\n";
        outfile << " [Cl- ] " << std::scientific << std::setprecision(3) << cl_print << "\n";
        outfile << " [NO3- ] " << std::scientific << std::setprecision(3) << no3_print << "\n";
        outfile << " [SO4-- ] " << std::scientific << std::setprecision(3) << so4_print << "\n";
        outfile << " [HSO4- ] " << std::scientific << std::setprecision(3) << hso4_print << "\n";
        outfile << " [Ca ] " << std::scientific << std::setprecision(3) << ca_print << "\n";
        outfile << " [K ] " << std::scientific << std::setprecision(3) << k_print << "\n";
        outfile << " [Mg ] " << std::scientific << std::setprecision(3) << mg_print << "\n";
        outfile << " [NH3 ] " << std::scientific << std::setprecision(3) << gnh3_print << "\n";
        outfile << " [HNO3 ] " << std::scientific << std::setprecision(3) << ghno3_print << "\n";
        outfile << " [HCL ] " << std::scientific << std::setprecision(3) << ghcl_print << "\n";
        outfile << " [Wat(NH4)2SO4] " << std::scientific << std::setprecision(3) << wat_nh42so4 << "\n";
        outfile << " [WatNH4NO3] " << std::scientific << std::setprecision(3) << wat_nh4no3 << "\n";
        outfile << " [WatOrg] " << std::scientific << std::setprecision(3) << wat_org << "\n";
        outfile << " [pH ] " << std::scientific << std::setprecision(3) << ph_print << "\n";
        outfile << " [IONIC STRENGTH] " << std::scientific << std::setprecision(3) << ionic_print << "\n";
        outfile << " ================================================\n\n";
    }

    std::cout << "Successfully generated " << out_file << std::endl;
    return 0;
}
