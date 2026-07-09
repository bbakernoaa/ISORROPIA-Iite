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

    while (std::getline(infile, line)) {
        std::string trimmed = trim_str(line);
        if (trimmed.empty() || trimmed[0] == 'C' || trimmed[0] == '#' || trimmed[0] == '*') continue;
        
        if (!headers_done) {
            // Locate species column header row
            if (trimmed.find("Na") != std::string::npos && trimmed.find("SO4") != std::string::npos) {
                headers_done = true;
            }
            continue;
        }

        std::stringstream ss(trimmed);
        InputRecord rec;
        if (ss >> rec.na >> rec.so4 >> rec.nh3 >> rec.no3 >> rec.cl >> rec.ca >> rec.k >> rec.mg >> rec.org >> rec.k_org >> rec.rho_org >> rec.rh >> rec.temp) {
            records.push_back(rec);
        }
    }

    std::cout << "Parsed " << records.size() << " records from input file." << std::endl;

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
        // Convert input masses from ug/m3 to mol/m3 inside the solver
        input.w[0] = std::max((rec.na / 23.0) * 1e-6, 0.0);
        input.w[1] = std::max((rec.so4 / 98.0) * 1e-6, 0.0);
        input.w[2] = std::max((rec.nh3 / 17.0) * 1e-6, 0.0);
        input.w[3] = std::max((rec.no3 / 63.0) * 1e-6, 0.0);
        input.w[4] = std::max((rec.cl / 36.5) * 1e-6, 0.0);
        input.w[5] = std::max((rec.ca / 40.1) * 1e-6, 0.0);
        input.w[6] = std::max((rec.k / 39.1) * 1e-6, 0.0);
        input.w[7] = std::max((rec.mg / 24.3) * 1e-6, 0.0);

        // Convert Organic concentrations to kilograms/m3 (ug/m3 to kg/m3) and density to kg/m3
        input.org[0] = rec.org * 1e-9;
        input.org[1] = rec.k_org;
        input.org[2] = rec.rho_org * 1e3;

        input.rh = rec.rh;
        input.temp = rec.temp;
        
        // Check if filename contains 'Reverse' or 'reverse'
        std::string lower_file = input_file;
        std::transform(lower_file.begin(), lower_file.end(), lower_file.begin(), ::tolower);
        if (lower_file.find("reverse") != std::string::npos) {
            input.iprob = 1; // Reverse Problem
        } else {
            input.iprob = 0; // Forward Problem
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
        
        double gnh3_print  = state.gnh3 * 17.0 * 1e6; // Gaseous Ammonia (MW=17) -> ug/m3
        double ghno3_print = state.ghno3 * 63.0 * 1e6; // Gaseous Nitric Acid (MW=63) -> ug/m3
        
        double ph_print    = (state.water > 1e-20 && state.molal[1] > 1e-30) ? -std::log10(state.molal[1] / state.water) : 7.0;
        double ionic_print = state.ionic;
        
        double wat_nh42so4 = state.watcmp[3] * 1e9;
        double wat_nh4no3  = state.watcmp[4] * 1e9;
        double wat_org     = state.watcmp[23] * 1e9;

        outfile << "Record " << r + 1 << " Solution:" << std::endl;
        outfile << " [WATER ] " << std::scientific << std::setprecision(3) << water_val << "\n";
        outfile << " [H+ ] " << std::scientific << std::setprecision(3) << h_print << "\n";
        outfile << " [NH4+ ] " << std::scientific << std::setprecision(3) << nh4_print << "\n";
        outfile << " [NO3- ] " << std::scientific << std::setprecision(3) << no3_print << "\n";
        outfile << " [SO4-- ] " << std::scientific << std::setprecision(3) << so4_print << "\n";
        outfile << " [HSO4- ] " << std::scientific << std::setprecision(3) << hso4_print << "\n";
        outfile << " [NH3 ] " << std::scientific << std::setprecision(3) << gnh3_print << "\n";
        outfile << " [HNO3 ] " << std::scientific << std::setprecision(3) << ghno3_print << "\n";
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
