import re
import os

def parse_fortran_zsr():
    """
    Parses uncommented DATA AW* statements from isolite1_0_src/isocom.f
    and generates C++ code with expanded values.
    """
    with open("isolite1_0_src/isocom.f", "r") as f:
        lines = f.readlines()

    # We want to find uncommented DATA blocks starting from line 500 to 870
    zsr_data = {}
    current_array = None
    accumulated_str = ""

    # Matches uncommented lines like:       DATA AWAS/10*187.72,
    # or continuation lines:      & 158.13,134.41,
    data_start_re = re.compile(r'^\s{6}DATA\s+(AW[A-Z0-9]+)\s*/', re.IGNORECASE)
    continuation_re = re.compile(r'^\s{5}([\&\$0-9])(.*)')

    for i, line in enumerate(lines):
        line_num = i + 1
        # Skip commented lines
        if line.startswith('C') or line.startswith('c') or line.startswith('*'):
            # If we were accumulating, the comment ends the data statement
            if current_array:
                zsr_data[current_array] = accumulated_str
                current_array = None
                accumulated_str = ""
            continue

        match_start = data_start_re.match(line)
        if match_start:
            # Save previous if any
            if current_array:
                zsr_data[current_array] = accumulated_str
            current_array = match_start.group(1).upper()
            # Everything after the slash
            parts = line.split('/', 1)
            accumulated_str = parts[1].strip()
            continue

        if current_array:
            match_cont = continuation_re.match(line)
            if match_cont:
                content = match_cont.group(2).strip()
                accumulated_str += " " + content
            else:
                # Line does not match continuation and is not a comment, so the DATA statement is done
                zsr_data[current_array] = accumulated_str
                current_array = None
                accumulated_str = ""

    # Ensure last one is saved
    if current_array:
        zsr_data[current_array] = accumulated_str

    # Now let's expand and format each array
    cpp_arrays = {}
    for name, raw_val_str in zsr_data.items():
        # Remove trailing slash and comments
        val_str = raw_val_str
        if '/' in val_str:
            val_str = val_str.split('/', 1)[0]
        val_str = val_str.replace('\n', ' ').strip()
        
        # Tokenize values by comma or space
        tokens = re.split(r'[\s,]+', val_str)
        expanded_values = []
        for token in tokens:
            if not token:
                continue
            # Check for repetition like 10*187.72 or 44*1000.00
            if '*' in token:
                rep_parts = token.split('*')
                count = int(rep_parts[0])
                val = rep_parts[1]
                # In Fortran, float representation like 1.e5 or 1.D5
                val = val.replace('d', 'e').replace('D', 'e')
                for _ in range(count):
                    expanded_values.append(float(val))
            else:
                val = token.replace('d', 'e').replace('D', 'e')
                expanded_values.append(float(val))

        if len(expanded_values) != 100:
            print(f"Warning: {name} has {len(expanded_values)} elements instead of 100!")
        cpp_arrays[name] = expanded_values

    return cpp_arrays

if __name__ == "__main__":
    arrays = parse_fortran_zsr()
    print(f"Parsed {len(arrays)} arrays successfully.")
    
    # Let's generate a C++ implementation function
    cpp_code = """#include "Isorropia/Solver.hpp"

namespace Isorropia {

void State::initialize_water_activities() {
"""
    for name, vals in sorted(arrays.items()):
        cpp_name = name.lower()
        cpp_code += f"    // {name} pure salt water activity grid (NZSR=100)\n"
        cpp_code += f"    {cpp_name} = {{\n"
        # Format values 10 per line
        for chunk_idx in range(0, 100, 10):
            chunk = vals[chunk_idx : chunk_idx + 10]
            chunk_str = ", ".join([f"{v:.6e}" for v in chunk])
            cpp_code += f"        {chunk_str}"
            if chunk_idx + 10 < 100:
                cpp_code += ",\n"
            else:
                cpp_code += "\n"
        cpp_code += "    };\n\n"

    # Add unused awcs array as zeroes explicitly or matching Fortran COMMON
    cpp_code += "    // AWCS is declared in COMMON /ZSR/ but unused in legacy code. Initializing to zero.\n"
    cpp_code += "    awcs.fill(0.0);\n"
    cpp_code += "}\n\n} // namespace Isorropia\n"

    with open("src/WaterActivities.cpp", "w") as f:
        f.write(cpp_code)
    print("Successfully generated src/WaterActivities.cpp")
