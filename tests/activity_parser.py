#!/usr/bin/env python3
import re


def parse_isocom(file_path):
    with open(file_path) as f:
        lines = f.readlines()

    blocks = {
        198: [],
        223: [],
        248: [],
        273: [],
        298: [],
        323: []
    }

    # Line ranges for BLOCK DATA KMCFxxx
    # Let's search for BLOCK DATA KMCFxxx and find the range of each block
    block_starts = []
    for i, line in enumerate(lines):
        m = re.search(r'BLOCK\s+DATA\s+KMCF(\d+)', line, re.IGNORECASE)
        if m:
            temp = int(m.group(1))
            block_starts.append((temp, i))

    block_starts.sort(key=lambda x: x[1])

    print("Found block starts:")
    for temp, idx in block_starts:
        print(f"  Temp: {temp}K starts at line {idx+1}")

    # For each block, extract arrays BNC01M to BNC23M
    # Each block goes from its start until the next block start, or END
    for idx_block, (temp, start_line_idx) in enumerate(block_starts):
        end_line_idx = len(lines)
        if idx_block + 1 < len(block_starts):
            end_line_idx = block_starts[idx_block + 1][1]

        block_lines = lines[start_line_idx:end_line_idx]
        
        # We want to extract DATA BNCxxM / ... /
        # We also want to capture comments immediately preceding DATA BNCxxM /
        current_comments = []
        i = 0
        while i < len(block_lines):
            line = block_lines[i].strip()
            
            # Check if it's a comment
            if line.startswith('C') or line.startswith('c') or line.startswith('*'):
                # Extract comment text
                comment_text = re.sub(r'^[Cc\*]\s*[\*]*\s*', '', block_lines[i]).strip()
                if comment_text:
                    current_comments.append(comment_text)
                i += 1
                continue
            
            # Check if it's a DATA BNCxxM line
            m = re.match(r'DATA\s+BNC(\d+)M\s*/', line, re.IGNORECASE)
            if m:
                arr_idx = int(m.group(1))
                # Now extract elements until '/'
                data_str = ""
                # Get the rest of this line after '/'
                rest_of_line = block_lines[i][block_lines[i].find('/') + 1:]
                data_str += rest_of_line
                
                # Read subsequent lines until we find '/'
                i += 1
                while i < len(block_lines) and '/' not in block_lines[i]:
                    line_to_add = block_lines[i]
                    # Strip continuation character '&' if present at start of line
                    line_to_add_stripped = line_to_add.strip()
                    if line_to_add_stripped.startswith('&'):
                        line_to_add_stripped = line_to_add_stripped[1:]
                    data_str += " " + line_to_add_stripped
                    i += 1
                
                if i < len(block_lines) and '/' in block_lines[i]:
                    # Add part before '/'
                    part_before = block_lines[i][:block_lines[i].find('/')]
                    line_to_add_stripped = part_before.strip()
                    if line_to_add_stripped.startswith('&'):
                        line_to_add_stripped = line_to_add_stripped[1:]
                    data_str += " " + line_to_add_stripped
                    i += 1
                
                # Split elements by comma and convert to float list
                # Some floats have spaces, newlines, etc.
                elements = []
                for val in data_str.replace('\n', ' ').split(','):
                    val = val.strip()
                    if val:
                        # Fortran floating point like .050, -.050 etc.
                        elements.append(float(val))
                
                species_name = " ".join(current_comments) if current_comments else f"Species {arr_idx}"
                blocks[temp].append({
                    'index': arr_idx,
                    'species': species_name,
                    'data': elements
                })
                current_comments = []
            else:
                if line != "":
                    # If it's some other non-empty line (not comment, not DATA), reset comments
                    # unless it's just continuation line which shouldn't happen outside DATA
                    current_comments = []
                i += 1

    return blocks

if __name__ == '__main__':
    isocom_path = "/Users/barry/Documents/ISORROPIA-Iite/isolite1_0_src/isocom.f"
    blocks = parse_isocom(isocom_path)
    
    # Check counts
    for temp, arrays in blocks.items():
        print(f"Temp {temp}K:")
        print(f"  Total arrays extracted: {len(arrays)}")
        for arr in arrays:
            if len(arr['data']) != 561:
                print(f"      ERROR: len is {len(arr['data'])} instead of 561!")

    # Generate src/ActivityCoefficients_KM.cpp
    print("\nGenerating src/ActivityCoefficients_KM.cpp...")
    cpp_code = """#include "Isorropia/Solver.hpp"
#include <cmath>
#include <algorithm>

namespace Isorropia {

"""
    # Write arrays for each temperature point
    for temp in sorted(blocks.keys()):
        arrays = blocks[temp]
        for arr in arrays:
            idx = arr['index']
            species = arr['species']
            cpp_code += f"// Temp {temp} K, Species BNC{idx:02d}M ({species})\n"
            cpp_code += f"static constexpr std::array<double, 561> bnc{idx:02d}m_{temp} = {{\n"
            # Format values 10 per line
            vals = arr['data']
            for chunk_idx in range(0, 561, 10):
                chunk = vals[chunk_idx : chunk_idx + 10]
                chunk_str = ", ".join([f"{v:.6e}" for v in chunk])
                cpp_code += f"    {chunk_str}"
                if chunk_idx + 10 < 561:
                    cpp_code += ",\n"
                else:
                    cpp_code += "\n"
            cpp_code += "};\n\n"

    # Write the main km_tab function
    cpp_code += """void State::km_tab(double ionic_strength, double temp_k, std::array<double, 23>& g0) {
    // 1. Calculate the nearest lookup position in 561-size grids
    int ipos = 0;
    if (ionic_strength <= 20.0) {
        ipos = std::min(static_cast<int>(std::round(20.0 * ionic_strength)) + 1, 400);
    } else {
        ipos = 400 + static_cast<int>(std::round(2.0 * ionic_strength - 40.0));
    }
    ipos = std::min(ipos, 561);
    size_t idx = static_cast<size_t>(ipos - 1);

    // 2. Determine temperature index IND (198 to 323 K in intervals of 25K)
    int ind = static_cast<int>((temp_k - 198.0) / 25.0 + 0.5) + 1;
    ind = std::max(1, std::min(ind, 6));

    // 3. Extract correct coefficients based on temperature index IND
    if (ind == 1) { // 198 K
"""
    # Map index 1 (198 K)
    for i in range(1, 24):
        cpp_code += f"        g0[{i-1}] = bnc{i:02d}m_198[idx];\n"
    
    # Map index 2 (223 K)
    cpp_code += "    } else if (ind == 2) { // 223 K\n"
    for i in range(1, 24):
        cpp_code += f"        g0[{i-1}] = bnc{i:02d}m_223[idx];\n"

    # Map index 3 (248 K)
    cpp_code += "    } else if (ind == 3) { // 248 K\n"
    for i in range(1, 24):
        cpp_code += f"        g0[{i-1}] = bnc{i:02d}m_248[idx];\n"

    # Map index 4 (273 K)
    cpp_code += "    } else if (ind == 4) { // 273 K\n"
    for i in range(1, 24):
        cpp_code += f"        g0[{i-1}] = bnc{i:02d}m_273[idx];\n"

    # Map index 5 (298 K)
    cpp_code += "    } else if (ind == 5) { // 298 K\n"
    for i in range(1, 24):
        cpp_code += f"        g0[{i-1}] = bnc{i:02d}m_298[idx];\n"

    # Map index 6 (323 K)
    cpp_code += "    } else { // 323 K\n"
    for i in range(1, 24):
        cpp_code += f"        g0[{i-1}] = bnc{i:02d}m_323[idx];\n"

    cpp_code += """    }
}

} // namespace Isorropia
"""

    with open("src/ActivityCoefficients_KM.cpp", "w") as f:
        f.write(cpp_code)
    print("Successfully generated src/ActivityCoefficients_KM.cpp")

